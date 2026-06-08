const api = {
  parameters: "/api/parameters",
  runs: "/api/runs",
  health: "/api/health",
};

const categoryOrder = [
  "Simulation",
  "Code",
  "MIMO",
  "Interleaving",
  "Decoder",
  "Architecture",
  "Channel",
  "General",
];

const quickNames = [
  "Min_SNR",
  "Max_SNR",
  "SNR_step",
  "errframe",
  "softiter",
  "Nh",
  "Nv",
  "Kh",
  "Kv",
  "Nt",
  "Nr",
  "trans_stream",
  "MIMODET",
  "polarconh",
  "polarconv",
  "MIMO_SVD_POWER_ALLOCATION_MODE",
];

const terminalStatuses = new Set(["completed", "failed", "terminated"]);

const state = {
  parameters: [],
  values: {},
  runs: [],
  selectedRunId: null,
  selectedRun: null,
  eventSource: null,
  metrics: new Map(),
  liveLineEl: null,
  snapshotMode: new URLSearchParams(window.location.search).has("snapshot"),
};

const els = {
  healthText: document.getElementById("healthText"),
  runName: document.getElementById("runName"),
  startBtn: document.getElementById("startBtn"),
  stopBtn: document.getElementById("stopBtn"),
  refreshBtn: document.getElementById("refreshBtn"),
  resetBtn: document.getElementById("resetBtn"),
  runCount: document.getElementById("runCount"),
  runsList: document.getElementById("runsList"),
  paramCount: document.getElementById("paramCount"),
  paramSearch: document.getElementById("paramSearch"),
  categoryFilter: document.getElementById("categoryFilter"),
  quickParams: document.getElementById("quickParams"),
  parameterGroups: document.getElementById("parameterGroups"),
  selectedRunTitle: document.getElementById("selectedRunTitle"),
  selectedRunMeta: document.getElementById("selectedRunMeta"),
  metricCount: document.getElementById("metricCount"),
  metricsBody: document.getElementById("metricsBody"),
  logView: document.getElementById("logView"),
  clearLogBtn: document.getElementById("clearLogBtn"),
  autoScroll: document.getElementById("autoScroll"),
};

document.addEventListener("DOMContentLoaded", async () => {
  wireEvents();
  await Promise.all([loadHealth(), loadParameters(), refreshRuns()]);
  setInterval(refreshRuns, 2500);
  drawIcons();
});

function wireEvents() {
  els.startBtn.addEventListener("click", startRun);
  els.stopBtn.addEventListener("click", () => terminateRun(state.selectedRunId));
  els.refreshBtn.addEventListener("click", async () => {
    await refreshRuns();
    await loadParameters(false);
  });
  els.resetBtn.addEventListener("click", () => loadParameters(true));
  els.paramSearch.addEventListener("input", renderParameters);
  els.categoryFilter.addEventListener("change", renderParameters);
  els.clearLogBtn.addEventListener("click", () => {
    els.logView.innerHTML = "";
    state.liveLineEl = null;
  });
}

async function loadHealth() {
  try {
    const data = await getJson(api.health);
    const msbuild = data.msbuild ? "MSBuild ready" : "MSBuild missing";
    els.healthText.textContent = `Release x64 · /O2 · ${msbuild}`;
  } catch {
    els.healthText.textContent = "Release x64 · /O2";
  }
}

async function loadParameters(resetValues = true) {
  const data = await getJson(api.parameters);
  state.parameters = data.parameters || [];
  if (resetValues) {
    state.values = {};
    for (const param of state.parameters) {
      state.values[param.name] = param.display_value;
    }
  } else {
    for (const param of state.parameters) {
      if (!(param.name in state.values)) {
        state.values[param.name] = param.display_value;
      }
    }
  }
  renderCategoryFilter();
  renderQuickParams();
  renderParameters();
  els.paramCount.textContent = `${state.parameters.length} active consts`;
  drawIcons();
}

function renderCategoryFilter() {
  const current = els.categoryFilter.value;
  const categories = [...new Set(state.parameters.map((param) => param.category))].sort(
    (a, b) => categoryRank(a) - categoryRank(b) || a.localeCompare(b),
  );
  els.categoryFilter.innerHTML = `<option value="">All</option>`;
  for (const category of categories) {
    const option = document.createElement("option");
    option.value = category;
    option.textContent = category;
    els.categoryFilter.appendChild(option);
  }
  if (categories.includes(current)) {
    els.categoryFilter.value = current;
  }
}

function renderQuickParams() {
  els.quickParams.innerHTML = "";
  const quick = quickNames
    .map((name) => state.parameters.find((param) => param.name === name))
    .filter(Boolean);
  for (const param of quick) {
    els.quickParams.appendChild(createParamField(param, true));
  }
}

function renderParameters() {
  const search = els.paramSearch.value.trim().toLowerCase();
  const category = els.categoryFilter.value;
  const groups = new Map();

  for (const param of state.parameters) {
    const haystack = `${param.name} ${param.ctype} ${param.comment}`.toLowerCase();
    if (search && !haystack.includes(search)) {
      continue;
    }
    if (category && param.category !== category) {
      continue;
    }
    if (!groups.has(param.category)) {
      groups.set(param.category, []);
    }
    groups.get(param.category).push(param);
  }

  els.parameterGroups.innerHTML = "";
  const entries = [...groups.entries()].sort(
    ([a], [b]) => categoryRank(a) - categoryRank(b) || a.localeCompare(b),
  );

  if (!entries.length) {
    const empty = document.createElement("div");
    empty.className = "empty-state";
    empty.textContent = "No parameters";
    els.parameterGroups.appendChild(empty);
    return;
  }

  for (const [group, params] of entries) {
    const section = document.createElement("section");
    section.className = "param-group";

    const title = document.createElement("div");
    title.className = "group-title";
    title.innerHTML = `<span>${escapeHtml(group)}</span><span class="group-count">${params.length}</span>`;
    section.appendChild(title);

    const grid = document.createElement("div");
    grid.className = "param-grid";
    for (const param of params) {
      grid.appendChild(createParamField(param, false));
    }
    section.appendChild(grid);
    els.parameterGroups.appendChild(section);
  }
}

function createParamField(param, quick) {
  const field = document.createElement("div");
  field.className = `param-field${quick ? " quick" : ""}${param.is_array && !quick ? " array" : ""}`;

  const label = document.createElement("div");
  label.className = "param-label";
  label.innerHTML = `
    <span class="param-name" title="${escapeHtml(param.name)}">${escapeHtml(param.name)}</span>
    <span class="param-type">${escapeHtml(typeLabel(param))}</span>
  `;
  field.appendChild(label);

  const control = createControl(param);
  field.appendChild(control);

  if (!quick && param.comment) {
    const comment = document.createElement("div");
    comment.className = "param-comment";
    comment.textContent = param.comment;
    field.appendChild(comment);
  }
  return field;
}

function createControl(param) {
  const value = state.values[param.name];
  if (param.kind === "bool") {
    const row = document.createElement("div");
    row.className = "toggle-row";
    const label = document.createElement("label");
    label.className = "toggle";
    const input = document.createElement("input");
    input.type = "checkbox";
    input.dataset.param = param.name;
    input.checked = Boolean(value);
    input.addEventListener("change", () => updateValue(param.name, input.checked, input));
    const span = document.createElement("span");
    label.append(input, span);
    row.appendChild(label);
    return row;
  }

  if (param.kind === "array") {
    const textarea = document.createElement("textarea");
    textarea.className = "param-textarea";
    textarea.dataset.param = param.name;
    textarea.value = value ?? "";
    textarea.spellcheck = false;
    textarea.addEventListener("input", () => updateValue(param.name, textarea.value, textarea));
    return textarea;
  }

  const input = document.createElement("input");
  input.className = "param-input";
  input.dataset.param = param.name;
  input.value = value ?? "";
  input.spellcheck = false;
  if (param.kind === "number") {
    input.type = "number";
    input.step = "any";
  } else {
    input.type = "text";
    if (param.choices?.length) {
      const listId = `choices-${param.name}-${Math.random().toString(16).slice(2)}`;
      const datalist = document.createElement("datalist");
      datalist.id = listId;
      for (const choice of param.choices) {
        const option = document.createElement("option");
        option.value = choice;
        datalist.appendChild(option);
      }
      input.setAttribute("list", listId);
      input.addEventListener("focus", () => {
        if (!document.getElementById(listId)) {
          document.body.appendChild(datalist);
        }
      });
      document.body.appendChild(datalist);
    }
  }
  input.addEventListener("input", () => updateValue(param.name, input.value, input));
  return input;
}

function updateValue(name, value, source) {
  state.values[name] = value;
  for (const control of document.querySelectorAll(`[data-param="${name}"]`)) {
    if (control === source) {
      continue;
    }
    if (control.type === "checkbox") {
      control.checked = Boolean(value);
    } else {
      control.value = value ?? "";
    }
  }
}

async function startRun() {
  els.startBtn.disabled = true;
  try {
    const payload = {
      name: els.runName.value.trim() || null,
      parameters: { ...state.values },
    };
    const run = await postJson(api.runs, payload);
    await refreshRuns();
    selectRun(run.id);
    els.runName.value = "";
  } catch (error) {
    appendUiLog(`Start failed: ${error.message}`, "stderr");
  } finally {
    els.startBtn.disabled = false;
  }
}

async function terminateRun(runId) {
  if (!runId) {
    return;
  }
  els.stopBtn.disabled = true;
  try {
    await postJson(`/api/runs/${runId}/terminate`, {});
    await refreshRuns();
  } catch (error) {
    appendUiLog(`Terminate failed: ${error.message}`, "stderr");
  }
}

async function refreshRuns() {
  const data = await getJson(api.runs);
  state.runs = data.runs || [];
  renderRuns();
  if (!state.selectedRunId && state.runs.length) {
    selectRun(state.runs[0].id);
    return;
  }
  if (state.selectedRunId) {
    const selected = state.runs.find((run) => run.id === state.selectedRunId);
    if (selected) {
      state.selectedRun = { ...(state.selectedRun || {}), ...selected };
      updateSelectedHeader(state.selectedRun);
    }
  }
}

function renderRuns() {
  els.runCount.textContent = state.runs.length;
  els.runsList.innerHTML = "";
  if (!state.runs.length) {
    const empty = document.createElement("div");
    empty.className = "empty-state";
    empty.textContent = "No runs";
    els.runsList.appendChild(empty);
    return;
  }

  for (const run of state.runs) {
    const card = document.createElement("button");
    card.className = `run-card${run.id === state.selectedRunId ? " selected" : ""}`;
    card.addEventListener("click", () => selectRun(run.id));

    const metricsCount = Array.isArray(run.metrics) ? run.metrics.length : 0;
    card.innerHTML = `
      <div class="run-title-row">
        <div>
          <div class="run-title">${escapeHtml(run.name)}</div>
          <div class="run-meta">${formatTime(run.created_at)} · ${metricsCount} pts</div>
        </div>
        ${statusPill(run.status)}
      </div>
      <div class="run-meta">${escapeHtml(run.message || "")}</div>
    `;

    const actions = document.createElement("div");
    actions.className = "run-actions";
    const pid = document.createElement("span");
    pid.className = "run-meta";
    pid.textContent = run.pid ? `PID ${run.pid}` : run.compile_pid ? `Build ${run.compile_pid}` : "";
    const stop = document.createElement("button");
    stop.className = "danger-btn small-btn";
    stop.title = "Terminate";
    stop.innerHTML = `<i data-lucide="square"></i>`;
    stop.disabled = terminalStatuses.has(run.status);
    stop.addEventListener("click", (event) => {
      event.stopPropagation();
      terminateRun(run.id);
    });
    actions.append(pid, stop);
    card.appendChild(actions);
    els.runsList.appendChild(card);
  }
  drawIcons();
}

async function selectRun(runId) {
  if (!runId || state.selectedRunId === runId) {
    return;
  }
  state.selectedRunId = runId;
  state.metrics = new Map();
  state.liveLineEl = null;
  els.logView.innerHTML = "";
  renderMetrics();
  renderRuns();

  const cachedRun = state.runs.find((run) => run.id === runId);
  if (cachedRun) {
    state.selectedRun = cachedRun;
    updateSelectedHeader(cachedRun);
    for (const metric of cachedRun.metrics || []) {
      state.metrics.set(metricKey(metric.ebn0), metric);
    }
    renderMetrics();
  }

  if (state.eventSource) {
    state.eventSource.close();
  }
  if (state.snapshotMode) {
    try {
      const run = await getJson(`/api/runs/${runId}`);
      state.selectedRun = run;
      hydrateRun(run);
    } catch (error) {
      appendUiLog(`Snapshot load failed: ${error.message}`, "stderr");
    }
    return;
  }
  const source = new EventSource(`/api/runs/${runId}/events`);
  state.eventSource = source;

  source.addEventListener("snapshot", (event) => {
    const run = JSON.parse(event.data);
    state.selectedRun = run;
    hydrateRun(run);
  });
  source.addEventListener("status", (event) => {
    const run = JSON.parse(event.data);
    state.selectedRun = { ...(state.selectedRun || {}), ...run };
    updateSelectedHeader(state.selectedRun);
    refreshRuns();
  });
  source.addEventListener("process", (event) => {
    const run = JSON.parse(event.data);
    state.selectedRun = { ...(state.selectedRun || {}), ...run };
    updateSelectedHeader(state.selectedRun);
  });
  source.addEventListener("log", (event) => {
    const item = JSON.parse(event.data);
    appendUiLog(item.text, item.stream);
  });
  source.addEventListener("log_replace", (event) => {
    const item = JSON.parse(event.data);
    replaceLiveLog(item.text, item.stream);
  });
  source.addEventListener("metric", (event) => {
    const item = JSON.parse(event.data);
    upsertMetric(item.metric);
  });
  source.onerror = () => {
    if (state.selectedRun) {
      els.selectedRunMeta.textContent = `${state.selectedRun.status} · event stream reconnecting`;
    }
  };
}

function hydrateRun(run) {
  updateSelectedHeader(run);
  els.logView.innerHTML = "";
  state.liveLineEl = null;
  for (const item of run.logs || []) {
    appendUiLog(item.text, item.stream, false);
  }
  if (run.live_line) {
    replaceLiveLog(run.live_line, "stdout");
  }
  state.metrics = new Map();
  for (const metric of run.metrics || []) {
    state.metrics.set(metricKey(metric.ebn0), metric);
  }
  renderMetrics();
  maybeScrollLog();
}

function updateSelectedHeader(run) {
  if (!run) {
    els.selectedRunTitle.textContent = "Output";
    els.selectedRunMeta.textContent = "No run selected";
    els.stopBtn.disabled = true;
    return;
  }
  els.selectedRunTitle.textContent = run.name || "Output";
  const bits = [run.status, run.message];
  if (run.pid) bits.push(`PID ${run.pid}`);
  els.selectedRunMeta.textContent = bits.filter(Boolean).join(" · ");
  els.stopBtn.disabled = terminalStatuses.has(run.status);
}

function appendUiLog(text, stream = "stdout", scroll = true) {
  if (state.liveLineEl) {
    state.liveLineEl.classList.remove("live");
    state.liveLineEl = null;
  }
  const line = document.createElement("div");
  line.className = `log-line ${stream || "stdout"}`;
  line.textContent = text ?? "";
  els.logView.appendChild(line);
  trimLogDom();
  if (scroll) maybeScrollLog();
}

function replaceLiveLog(text, stream = "stdout") {
  if (!state.liveLineEl) {
    state.liveLineEl = document.createElement("div");
    state.liveLineEl.className = `log-line live ${stream || "stdout"}`;
    els.logView.appendChild(state.liveLineEl);
  }
  state.liveLineEl.textContent = text ?? "";
  maybeScrollLog();
}

function trimLogDom() {
  while (els.logView.children.length > 2200) {
    els.logView.removeChild(els.logView.firstChild);
  }
}

function maybeScrollLog() {
  if (els.autoScroll.checked) {
    els.logView.scrollTop = els.logView.scrollHeight;
  }
}

function upsertMetric(metric) {
  state.metrics.set(metricKey(metric.ebn0), metric);
  renderMetrics();
}

function renderMetrics() {
  const metrics = [...state.metrics.values()].sort((a, b) => Number(a.ebn0) - Number(b.ebn0));
  els.metricCount.textContent = metrics.length;
  els.metricsBody.innerHTML = "";
  if (!metrics.length) {
    els.metricsBody.innerHTML = `<tr class="empty-row"><td colspan="5">Waiting for data</td></tr>`;
    return;
  }
  for (const metric of metrics) {
    const tr = document.createElement("tr");
    tr.innerHTML = `
      <td>${formatNumber(metric.ebn0, 2)}</td>
      <td>${formatRate(metric.fer)}</td>
      <td>${formatRate(metric.ber)}</td>
      <td>${formatCount(metric.frames)}</td>
      <td>${formatCount(metric.bits)}</td>
    `;
    els.metricsBody.appendChild(tr);
  }
}

async function getJson(url) {
  const response = await fetch(url);
  if (!response.ok) {
    throw new Error(await response.text());
  }
  return response.json();
}

async function postJson(url, payload) {
  const response = await fetch(url, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  });
  if (!response.ok) {
    throw new Error(await response.text());
  }
  return response.json();
}

function typeLabel(param) {
  return param.is_array ? `${param.ctype}${param.array_suffix}` : param.ctype;
}

function categoryRank(category) {
  const index = categoryOrder.indexOf(category);
  return index >= 0 ? index : 99;
}

function statusPill(status) {
  return `<span class="status-pill status-${escapeHtml(status)}">${escapeHtml(status)}</span>`;
}

function metricKey(value) {
  return Number(value).toFixed(6);
}

function formatRate(value) {
  if (value === undefined || value === null) {
    return "-";
  }
  const number = Number(value);
  if (!Number.isFinite(number)) {
    return "-";
  }
  if (number > 0 && number < 0.001) {
    return number.toExponential(3);
  }
  return number.toFixed(6);
}

function formatNumber(value, digits = 2) {
  const number = Number(value);
  if (!Number.isFinite(number)) {
    return "-";
  }
  return number.toFixed(digits);
}

function formatCount(value) {
  if (value === undefined || value === null) {
    return "-";
  }
  const number = Number(value);
  if (!Number.isFinite(number)) {
    return "-";
  }
  return number.toLocaleString();
}

function formatTime(seconds) {
  if (!seconds) {
    return "";
  }
  return new Date(seconds * 1000).toLocaleTimeString();
}

function escapeHtml(value) {
  return String(value ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;");
}

function drawIcons() {
  if (window.lucide) {
    window.lucide.createIcons();
  }
}
