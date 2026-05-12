# include"mat_io.h"
# include<Eigen/Dense>
# include<string>
# include<iomanip>
void save_mat_txt(const Eigen::MatrixXf& mat, const std::string& file_name) {
    std::ofstream file(file_name);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << file_name << std::endl;
        return;
    }

    // 获取矩阵的行数和列数
    int rows = mat.rows();
    int cols = mat.cols();

    // 将矩阵的数据写入文件
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            file << mat(i, j);
            if (j < cols - 1) {
                file << " "; // 使用空格分隔
            }
        }
        file << std::endl; // 换行
    }

    file.close(); // 关闭文件
}

void save_mat_txt(const float* mat, int cols, const std::string& file_name) {
    std::ofstream file(file_name);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << file_name << std::endl;
        return;
    }

 
    // 将矩阵的数据写入文件
    for (int i = 0; i < cols; i++)
    {
        file << mat[i];
        file << "  ";
    }
    file << endl;
    file.close(); // 关闭文件
}

void save_array_txt(int* arr, int len, const std::string& file_name) {
    std::ofstream ofile(file_name);
    for (int i = 0; i < len; i++)
    {
        ofile << arr[i] << "  ";
    }
}

void save_array_txt_commasep(int* arr, int len, const std::string& file_name) {
    std::ofstream ofile(file_name);
    for (int i = 0; i < len; i++)
    {
        ofile << arr[i] << ",";
    }
}

void save_array_txt(float SNR, float BER, float FER, const std::string& file_name) {
    std::ofstream file(file_name, std::ios::app);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << file_name << std::endl;
        return;
    }

    file << setw(20) << SNR << setw(20) << BER << setw(20) << FER << "\n";

    file.close(); // 关闭文件
}

void save_array_txt(float SNR, float FER, const std::string& file_name) {
    std::ofstream file(file_name, std::ios::app);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << file_name << std::endl;
        return;
    }

    file << setw(20) << SNR << setw(20)  << FER << "\n";

    file.close(); // 关闭文件
}

void display_array(int** arr, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            std::cout << arr[i][j] << " ";
        }
        std::cout << std::endl;
    }
}


