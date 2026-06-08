from __future__ import annotations

import asyncio
import codecs
import html
import json
import locale
import os
import re
import shutil
import subprocess
import time
import uuid
from collections import deque
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel


UI_DIR = Path(__file__).resolve().parent
ROOT_DIR = UI_DIR.parent
PROJECT_DIR = ROOT_DIR / "platform_1"
SETTING_PATH = PROJECT_DIR / "setting.h"
SOLUTION_PATH = ROOT_DIR / "platform_1.sln"
BUILD_EXE_PATH = ROOT_DIR / "x64" / "Release" / "platform_1.exe"
RUNS_DIR = UI_DIR / "runs"

CREATE_NO_WINDOW = getattr(subprocess, "CREATE_NO_WINDOW", 0)
CREATE_NEW_PROCESS_GROUP = getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)

RUNTIME_SKIP_DIRS = {".vs", "Debug", "Release", "x64"}
RUNTIME_OUTPUT_PREFIXES = (
    "AllPowerAllocResult",
    "AllRateRegularResult",
    "EbN0VsCodeRate",
    "OptimizationResult",
    "all_errors_log",
    "crc_check",
    "crc_replacement_codeword_log",
    "dec_timing",
    "det_timing",
    "error_pattern_elaa",
    "error_patterns",
    "index_output",
    "progress_iter",
    "result",
    "run_err",
    "run_log",
    "sysresult",
    "temp_output",
    "total_frames",
    "total_timing",
    "y_out",
)

DECL_RE = re.compile(
    r"^(?P<indent>\s*)const\s+"
    r"(?P<ctype>(?:std::)?string|int|float|double|bool)\s+"
    r"(?P<name>[A-Za-z_]\w*)\s*"
    r"(?P<array>(?:\[[^\]]+\])*)\s*=\s*"
    r"(?P<rest>.+)$"
)

LIVE_METRIC_RE = re.compile(
    r"SNR\s*=\s*(?P<snr>[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"
    r".*?FER\s*=\s*(?P<fer_err>\d+)\s*/\s*(?P<fer_total>\d+)\s*=\s*"
    r"(?P<fer>[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"
    r".*?BER\s*=\s*(?P<ber_err>\d+)\s*/\s*(?P<ber_total>\d+)\s*=\s*"
    r"(?P<ber>[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)",
    re.IGNORECASE,
)
SUMMARY_ROW_RE = re.compile(
    r"^\s*(?P<snr>[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s+"
    r"(?P<value>[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*$"
)

CHOICES: dict[str, list[str]] = {
    "polarcon": ["RM", "5G", "beta", "GA"],
    "polarconh": ["RM-GA", "RM", "5G", "beta", "GA"],
    "polarconv": ["RM", "5G", "beta", "GA"],
    "MOD_DIM": ["Temporal", "Spatial"],
    "MIMODET": ["MMSE", "KBEST", "EP"],
    "softmethod": ["Pyndiah", "MITSO"],
    "polarde": ["SLD", "SCL"],
    "system_archi": ["SDD", "JDD"],
    "MIMO_SVD_POWER_ALLOCATION_MODE": ["SEPARATE", "UNIFORM", "WATERFILLING"],
}


class StartRunRequest(BaseModel):
    name: str | None = None
    parameters: dict[str, Any] = {}


class Parameter(BaseModel):
    name: str
    ctype: str
    value: str
    display_value: Any
    kind: str
    is_array: bool
    array_suffix: str
    category: str
    comment: str
    line: int
    choices: list[str] = []


class RunState:
    def __init__(self, run_id: str, name: str, parameters: dict[str, Any]) -> None:
        self.id = run_id
        self.name = name
        self.parameters = parameters
        self.status = "queued"
        self.created_at = time.time()
        self.updated_at = self.created_at
        self.finished_at: float | None = None
        self.message = "Queued"
        self.logs: deque[dict[str, Any]] = deque(maxlen=6000)
        self.events: list[dict[str, Any]] = []
        self.live_line = ""
        self.metrics: dict[str, dict[str, Any]] = {}
        self.summary_mode: str | None = None
        self.process: asyncio.subprocess.Process | None = None
        self.compile_process: asyncio.subprocess.Process | None = None
        self.terminate_requested = False
        self.run_dir = RUNS_DIR / run_id
        self.work_dir = self.run_dir / "work"
        self.exe_path = self.run_dir / "platform_1.exe"

    def to_public(self, include_logs: bool = False) -> dict[str, Any]:
        data = {
            "id": self.id,
            "name": self.name,
            "status": self.status,
            "message": self.message,
            "created_at": self.created_at,
            "updated_at": self.updated_at,
            "finished_at": self.finished_at,
            "pid": self.process.pid if self.process else None,
            "compile_pid": self.compile_process.pid if self.compile_process else None,
            "live_line": self.live_line,
            "metrics": sorted(self.metrics.values(), key=lambda item: item["ebn0"]),
            "parameters": self.parameters,
        }
        if include_logs:
            data["logs"] = list(self.logs)
        return data

    def emit(self, event_type: str, payload: dict[str, Any]) -> None:
        self.updated_at = time.time()
        event = {
            "type": event_type,
            "time": self.updated_at,
            **payload,
        }
        self.events.append(event)

    def set_status(self, status: str, message: str) -> None:
        self.status = status
        self.message = message
        if status in {"completed", "failed", "terminated"}:
            self.finished_at = time.time()
        self.emit("status", self.to_public(include_logs=False))

    def add_log(self, text: str, stream: str = "stdout", replace: bool = False) -> None:
        clean = text.rstrip("\n")
        if replace:
            self.live_line = clean
            self.emit("log_replace", {"text": clean, "stream": stream})
        else:
            if clean:
                item = {"time": time.time(), "text": clean, "stream": stream}
                self.logs.append(item)
                self.emit("log", item)
            self.live_line = ""
        self.consume_metric(clean)

    def consume_metric(self, line: str) -> None:
        text = line.strip()
        if not text:
            return

        live = LIVE_METRIC_RE.search(text)
        if live:
            metric = {
                "ebn0": float(live.group("snr")),
                "fer": float(live.group("fer")),
                "ber": float(live.group("ber")),
                "frame_errors": int(live.group("fer_err")),
                "frames": int(live.group("fer_total")),
                "bit_errors": int(live.group("ber_err")),
                "bits": int(live.group("ber_total")),
            }
            self.upsert_metric(metric)
            return

        header = re.sub(r"\s+", " ", text).upper()
        if header.startswith("SNR FER"):
            self.summary_mode = "fer"
            return
        if header.startswith("SNR BER"):
            self.summary_mode = "ber"
            return

        if self.summary_mode:
            row = SUMMARY_ROW_RE.match(text)
            if row:
                ebn0 = float(row.group("snr"))
                value = float(row.group("value"))
                current = self.metrics.get(metric_key(ebn0), {"ebn0": ebn0})
                current[self.summary_mode] = value
                self.upsert_metric(current)

    def upsert_metric(self, metric: dict[str, Any]) -> None:
        key = metric_key(float(metric["ebn0"]))
        existing = self.metrics.get(key, {"ebn0": float(metric["ebn0"])})
        existing.update(metric)
        self.metrics[key] = existing
        self.emit("metric", {"metric": existing})


app = FastAPI(title="MIMO 2D Simulation UI")
app.mount("/static", StaticFiles(directory=UI_DIR / "static"), name="static")

runs: dict[str, RunState] = {}
compile_lock = asyncio.Lock()


@app.get("/")
async def index() -> FileResponse:
    return FileResponse(UI_DIR / "static" / "index.html")


@app.get("/api/health")
async def health() -> dict[str, Any]:
    return {
        "root": str(ROOT_DIR),
        "setting": str(SETTING_PATH),
        "solution": str(SOLUTION_PATH),
        "msbuild": find_msbuild(),
    }


@app.get("/api/parameters")
async def get_parameters() -> dict[str, Any]:
    params = parse_setting_parameters()
    return {"parameters": [param.model_dump() for param in params]}


@app.get("/api/runs")
async def list_runs() -> dict[str, Any]:
    ordered = sorted(runs.values(), key=lambda run: run.created_at, reverse=True)
    return {"runs": [run.to_public(include_logs=False) for run in ordered]}


@app.post("/api/runs")
async def start_run(request: StartRunRequest) -> dict[str, Any]:
    run_id = time.strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:8]
    name = request.name or f"Run {run_id[-8:]}"
    run = RunState(run_id, name, request.parameters)
    runs[run_id] = run
    run.run_dir.mkdir(parents=True, exist_ok=True)
    asyncio.create_task(run_lifecycle(run))
    return run.to_public(include_logs=True)


@app.get("/api/runs/{run_id}")
async def get_run(run_id: str) -> dict[str, Any]:
    run = runs.get(run_id)
    if not run:
        raise HTTPException(status_code=404, detail="Run not found")
    return run.to_public(include_logs=True)


@app.post("/api/runs/{run_id}/terminate")
async def terminate_run(run_id: str) -> dict[str, Any]:
    run = runs.get(run_id)
    if not run:
        raise HTTPException(status_code=404, detail="Run not found")
    run.terminate_requested = True
    run.add_log("Termination requested by UI.", stream="system")

    if run.process and run.process.returncode is None:
        await kill_process_tree(run.process.pid)
    if run.compile_process and run.compile_process.returncode is None:
        await kill_process_tree(run.compile_process.pid)
    if run.status in {"queued", "compiling", "running"}:
        run.set_status("terminated", "Terminated")
    return run.to_public(include_logs=False)


@app.get("/api/runs/{run_id}/events")
async def run_events(run_id: str) -> StreamingResponse:
    run = runs.get(run_id)
    if not run:
        raise HTTPException(status_code=404, detail="Run not found")

    async def event_stream():
        yield to_sse("snapshot", run.to_public(include_logs=True))
        index = len(run.events)
        last_heartbeat = time.time()
        while True:
            if index < len(run.events):
                for event in run.events[index:]:
                    yield to_sse(event["type"], event)
                index = len(run.events)
                last_heartbeat = time.time()
            elif time.time() - last_heartbeat > 12:
                yield ": heartbeat\n\n"
                last_heartbeat = time.time()
            await asyncio.sleep(0.25)

    return StreamingResponse(event_stream(), media_type="text/event-stream")


async def run_lifecycle(run: RunState) -> None:
    try:
        async with compile_lock:
            if run.terminate_requested:
                run.set_status("terminated", "Terminated before compile")
                return
            run.set_status("compiling", "Applying settings and compiling Release x64")
            run.add_log("Applying selected parameters to platform_1\\setting.h", stream="system")
            apply_parameter_overrides(run.parameters, run.run_dir)

            msbuild = find_msbuild()
            if not msbuild:
                run.set_status("failed", "MSBuild.exe was not found")
                return

            command = [
                msbuild,
                str(SOLUTION_PATH),
                "/m",
                "/p:Configuration=Release",
                "/p:Platform=x64",
                "/p:PreferredToolArchitecture=x64",
            ]
            run.add_log(" ".join(command), stream="system")
            compile_code = await run_subprocess(command, ROOT_DIR, run, phase="compile")
            run.compile_process = None
            if run.terminate_requested:
                run.set_status("terminated", "Terminated during compile")
                return
            if compile_code != 0:
                run.set_status("failed", f"Compile failed with exit code {compile_code}")
                return
            if not BUILD_EXE_PATH.exists():
                run.set_status("failed", f"Expected executable not found: {BUILD_EXE_PATH}")
                return
            shutil.copy2(BUILD_EXE_PATH, run.exe_path)
            run.add_log(f"Executable snapshot: {run.exe_path}", stream="system")
            prepare_run_workdir(run)
            run.add_log(f"Isolated working directory: {run.work_dir}", stream="system")

        if run.terminate_requested:
            run.set_status("terminated", "Terminated before run")
            return
        run.set_status("running", "Running simulation")
        code = await run_subprocess([str(run.exe_path)], run.work_dir, run, phase="run")
        run.process = None
        if run.terminate_requested:
            run.set_status("terminated", "Terminated")
        elif code == 0:
            run.set_status("completed", "Completed")
        else:
            run.set_status("failed", f"Process exited with code {code}")
    except Exception as exc:
        run.add_log(f"{type(exc).__name__}: {exc}", stream="stderr")
        run.set_status("failed", "Unexpected UI backend error")


async def run_subprocess(
    command: list[str],
    cwd: Path,
    run: RunState,
    phase: str,
) -> int:
    process = await asyncio.create_subprocess_exec(
        *command,
        cwd=str(cwd),
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        creationflags=CREATE_NO_WINDOW | CREATE_NEW_PROCESS_GROUP,
    )
    if phase == "compile":
        run.compile_process = process
    else:
        run.process = process
    run.emit("process", run.to_public(include_logs=False))

    await asyncio.gather(
        stream_reader(process.stdout, run, "stdout"),
        stream_reader(process.stderr, run, "stderr"),
    )
    return await process.wait()


async def stream_reader(
    stream: asyncio.StreamReader | None,
    run: RunState,
    stream_name: str,
) -> None:
    if stream is None:
        return
    encoding = locale.getpreferredencoding(False) or "utf-8"
    decoder = codecs.getincrementaldecoder(encoding)(errors="replace")
    buffer = ""
    while True:
        chunk = await stream.read(512)
        if not chunk:
            break
        text = decoder.decode(chunk)
        for char in text:
            if char == "\r":
                if buffer:
                    run.add_log(buffer, stream=stream_name, replace=True)
                    buffer = ""
            elif char == "\n":
                run.add_log(buffer, stream=stream_name, replace=False)
                buffer = ""
            else:
                buffer += char
    tail = decoder.decode(b"", final=True)
    if tail:
        buffer += tail
    if buffer:
        run.add_log(buffer, stream=stream_name, replace=False)


async def kill_process_tree(pid: int) -> None:
    if os.name == "nt":
        killer = await asyncio.create_subprocess_exec(
            "taskkill",
            "/PID",
            str(pid),
            "/T",
            "/F",
            stdout=asyncio.subprocess.DEVNULL,
            stderr=asyncio.subprocess.DEVNULL,
            creationflags=CREATE_NO_WINDOW,
        )
        await killer.wait()
    else:
        try:
            os.kill(pid, 15)
        except ProcessLookupError:
            pass


def parse_setting_parameters() -> list[Parameter]:
    lines = SETTING_PATH.read_text(encoding="utf-8", errors="replace").splitlines()
    params: list[Parameter] = []
    current_category = "General"
    for idx, line in enumerate(lines, start=1):
        stripped = line.strip()
        if stripped.startswith("//"):
            current_category = category_from_comment(stripped, current_category)
            continue
        parsed = parse_const_line(line)
        if not parsed:
            continue
        kind = infer_kind(parsed["ctype"], parsed["value"], bool(parsed["array"]))
        display_value = display_for_value(parsed["ctype"], parsed["value"], kind)
        param = Parameter(
            name=parsed["name"],
            ctype=parsed["ctype"],
            value=parsed["value"],
            display_value=display_value,
            kind=kind,
            is_array=bool(parsed["array"]),
            array_suffix=parsed["array"],
            category=categorize(parsed["name"], current_category),
            comment=parsed["comment"],
            line=idx,
            choices=CHOICES.get(parsed["name"], []),
        )
        params.append(param)
    return params


def parse_const_line(line: str) -> dict[str, str] | None:
    match = DECL_RE.match(line)
    if not match:
        return None
    value_part, comment = split_cpp_comment(match.group("rest"))
    if ";" not in value_part:
        return None
    value = value_part.rsplit(";", 1)[0].strip()
    return {
        "indent": match.group("indent"),
        "ctype": match.group("ctype"),
        "name": match.group("name"),
        "array": match.group("array") or "",
        "value": value,
        "comment": comment.strip(),
    }


def split_cpp_comment(text: str) -> tuple[str, str]:
    in_string = False
    escaped = False
    for index in range(len(text) - 1):
        char = text[index]
        if char == "\\" and in_string:
            escaped = not escaped
            continue
        if char == '"' and not escaped:
            in_string = not in_string
        escaped = False
        if not in_string and text[index : index + 2] == "//":
            return text[:index].rstrip(), text[index + 2 :].strip()
    return text.rstrip(), ""


def category_from_comment(comment: str, current: str) -> str:
    text = comment.strip("/ ").lower()
    if "mimo" in text:
        return "MIMO"
    if "soft" in text or "iteration" in text:
        return "Decoder"
    if "simulation" in text:
        return "Simulation"
    if "siso" in text or "sdd" in text or "jdd" in text:
        return "Architecture"
    if "coding" in text or "1d" in text or "2d" in text:
        return "Code"
    if "other" in text:
        return "Decoder"
    return current


def categorize(name: str, fallback: str) -> str:
    simulation = {"Min_SNR", "Max_SNR", "SNR_step", "errframe", "record_frames", "weirditer"}
    code = {
        "dimensions",
        "N",
        "K",
        "Nh",
        "Nv",
        "Kh",
        "Kv",
        "Nv_1D",
        "polarcon",
        "polarconh",
        "polarconv",
        "sys",
        "irregular_coding",
        "num_time_constructions",
        "K_time_array",
        "border",
        "amplify_ratio",
        "amplify_ratio_array",
        "irregular_iteration",
    }
    interleaving = {
        "flag_interleave",
        "flag_interleave_all",
        "flag_interleave_partial",
        "flag_interleave_spatial",
        "flag_interleave_block",
        "flag_interleave_inner",
        "interleave_rows",
        "flag_separate_interleaving",
    }
    mimo_prefixes = ("mimo_", "MIMO_", "Nt", "Nr", "trans_stream", "ModType", "MOD_", "MIMODET", "K_KBEST", "iter_MIMO")
    decoder = {
        "softiter",
        "softmethod",
        "SOalpha",
        "alphalen",
        "alpha",
        "usebeta",
        "adaptive_beta",
        "betalen",
        "beta",
        "a_beta",
        "b_beta",
        "k_beta",
        "blk_cnt",
        "CRCL",
        "CRC_CHECK_ENABLE",
        "PARITY_CHECK_ENABLE",
        "CRC_VERTICAL",
        "CRC_HORIZONTAL",
        "L",
        "L_WHOLE_CODEWORD",
        "L_FINAL_SEARCH",
        "L_SUB_CODEWORD",
        "L_MIMO_BLK",
        "polarde",
        "polarde_array",
        "half_iter",
        "bit_flip_pyndiah",
        "double_stage",
        "L_LSD",
        "p",
        "horizontal_ml",
        "vertical_ml",
        "max_llr_ratio",
        "L_max",
        "disable_max_ratio",
        "test_sc",
        "test_cadsl",
    }
    architecture = {"atten", "system_archi", "beta_JDD", "restrict_llr", "showllr", "showllr_irr"}
    channel = {
        "non_uniform_channel",
        "USE_PLATFORM_CHANNEL",
        "USE_ELAA_CHANNEL",
        "N_FREQUENCY_SAMPLES",
        "diff_gain",
        "predefined_non_uniform",
        "ENUMERATE_CODEWORD",
        "channel_index_predefined",
        "EXCLUSION_NUM",
        "exclusion_set",
        "PROTECTED_PARITY_NUM",
        "protected_parity_set",
        "SPECIFIC_EXCLUSION",
        "OUTPUT_ERROR_PATTERN",
        "CLOG_RESTRICTION",
        "flag_information_mapping",
    }

    if name in simulation:
        return "Simulation"
    if name in code:
        return "Code"
    if name in interleaving:
        return "Interleaving"
    if name in decoder:
        return "Decoder"
    if name in architecture:
        return "Architecture"
    if name in channel:
        return "Channel"
    if name.startswith(mimo_prefixes):
        return "MIMO"
    return fallback or "General"


def infer_kind(ctype: str, value: str, is_array: bool) -> str:
    if is_array:
        return "array"
    low = value.strip().lower()
    if ctype == "bool" or low in {"true", "false"}:
        return "bool"
    if ctype in {"int", "float", "double"} and re.fullmatch(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", value.strip()):
        return "number"
    if ctype in {"string", "std::string"}:
        return "string"
    return "expression"


def display_for_value(ctype: str, value: str, kind: str) -> Any:
    if kind == "bool":
        return value.strip().lower() in {"true", "1"}
    if kind == "number":
        try:
            return float(value) if "." in value or "e" in value.lower() else int(value)
        except ValueError:
            return value
    if kind == "string":
        return unquote_cpp_string(value)
    if kind == "array":
        return strip_outer_braces(value)
    return value


def apply_parameter_overrides(overrides: dict[str, Any], run_dir: Path) -> None:
    original = SETTING_PATH.read_text(encoding="utf-8", errors="replace")
    (run_dir / "setting.before.h").write_text(original, encoding="utf-8")
    lines = original.splitlines(keepends=True)
    param_map = {param.name: param for param in parse_setting_parameters()}
    rendered: list[str] = []

    for line in lines:
        ending = "\n" if line.endswith("\n") else ""
        body = line[:-1] if ending else line
        parsed = parse_const_line(body)
        if parsed and parsed["name"] in overrides and parsed["name"] in param_map:
            param = param_map[parsed["name"]]
            formatted = format_cpp_value(param, overrides[parsed["name"]])
            comment = f" // {parsed['comment']}" if parsed["comment"] else ""
            body = f"{parsed['indent']}const {parsed['ctype']} {parsed['name']}{parsed['array']} = {formatted};{comment}"
        rendered.append(body + ending)

    generated = "".join(rendered)
    (run_dir / "setting.generated.h").write_text(generated, encoding="utf-8")
    SETTING_PATH.write_text(generated, encoding="utf-8")


def prepare_run_workdir(run: RunState) -> None:
    if run.work_dir.exists():
        shutil.rmtree(run.work_dir)
    run.work_dir.mkdir(parents=True, exist_ok=True)

    for directory in PROJECT_DIR.rglob("*"):
        if not directory.is_dir():
            continue
        relative = directory.relative_to(PROJECT_DIR)
        if relative.parts and relative.parts[0] in RUNTIME_SKIP_DIRS:
            continue
        (run.work_dir / relative).mkdir(parents=True, exist_ok=True)

    for file_path in PROJECT_DIR.iterdir():
        if not file_path.is_file():
            continue
        destination = run.work_dir / file_path.name
        if file_path.name == "setting.h":
            shutil.copy2(run.run_dir / "setting.generated.h", destination)
            continue
        if is_runtime_output_file(file_path.name):
            continue
        try:
            os.link(file_path, destination)
        except OSError:
            shutil.copy2(file_path, destination)


def is_runtime_output_file(name: str) -> bool:
    stem = Path(name).stem
    return any(stem.startswith(prefix) for prefix in RUNTIME_OUTPUT_PREFIXES)


def format_cpp_value(param: Parameter, value: Any) -> str:
    if value is None:
        return param.value
    if param.kind == "bool":
        if isinstance(value, bool):
            return "true" if value else "false"
        text = str(value).strip().lower()
        if text in {"1", "true", "yes", "on"}:
            return "true"
        if text in {"0", "false", "no", "off"}:
            return "false"
        return str(value).strip()
    if param.kind == "string":
        text = str(value).strip()
        if text.startswith('"') and text.endswith('"'):
            return text
        return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'
    if param.kind == "array":
        text = str(value).strip()
        if text.startswith("{") and text.endswith("}"):
            return text
        return "{ " + text + " }"
    return str(value).strip()


def strip_outer_braces(value: str) -> str:
    text = value.strip()
    if text.startswith("{") and text.endswith("}"):
        return text[1:-1].strip()
    return text


def unquote_cpp_string(value: str) -> str:
    text = value.strip()
    if len(text) >= 2 and text[0] == '"' and text[-1] == '"':
        return bytes(text[1:-1], "utf-8").decode("unicode_escape")
    return text


def metric_key(value: float) -> str:
    return f"{value:.8g}"


def find_msbuild() -> str | None:
    found = shutil.which("msbuild")
    if found:
        return found
    candidates = [
        Path(r"C:\Program Files\Microsoft Visual Studio\2022\Community\MSBuild\Current\Bin\amd64\MSBuild.exe"),
        Path(r"C:\Program Files\Microsoft Visual Studio\2022\Professional\MSBuild\Current\Bin\amd64\MSBuild.exe"),
        Path(r"C:\Program Files\Microsoft Visual Studio\2022\Enterprise\MSBuild\Current\Bin\amd64\MSBuild.exe"),
        Path(r"C:\Program Files\Microsoft Visual Studio\2022\BuildTools\MSBuild\Current\Bin\amd64\MSBuild.exe"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    return None


def to_sse(event: str, data: dict[str, Any]) -> str:
    return f"event: {event}\ndata: {json.dumps(data, ensure_ascii=False)}\n\n"


@app.get("/api/debug/setting-preview", include_in_schema=False)
async def setting_preview() -> dict[str, Any]:
    params = parse_setting_parameters()
    return {"count": len(params), "names": [param.name for param in params]}


def escape_for_log(text: str) -> str:
    return html.escape(text, quote=False)
