#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>
static constexpr std::double_t Sp =                    0.004;
static constexpr std::double_t Sd1 =                   0.01;
static constexpr std::double_t Sd2 =                   0.015;
static constexpr std::double_t Sd3 =                   0.015;
static constexpr std::double_t Sdk =                   0.025;
static constexpr std::double_t p_to_cpu =              0.15;
static constexpr std::double_t p_to_1disk =            0.2;
static constexpr std::double_t p_to_2disk =            0.15;
static constexpr std::double_t p_to_3disk =            0.1;
static constexpr std::double_t p_to_kdisks =           0.4;
static constexpr std::double_t p_disk_to_cpu =         0.35;
static constexpr std::double_t p_disk_to_same_disk =   0.25;
static constexpr std::double_t p_disk_to_k_disk =      0.4;
static constexpr std::uint64_t arr_of_k_vals[4] = {2, 3, 4, 5};
static constexpr std::double_t arr_of_r_vals[4] = {0.3, 0.55, 0.8, 1.0};
static constexpr std::double_t X1_const = 1 / ((1-(p_disk_to_cpu/(1-p_to_cpu))*((p_to_1disk+p_to_2disk+p_to_3disk)/(1-p_disk_to_same_disk)))*(1-p_to_cpu));
static constexpr std::double_t X2_const = X1_const * (p_to_1disk/(1-p_disk_to_same_disk));
static constexpr std::double_t X3_const = X1_const * (p_to_2disk/(1-p_disk_to_same_disk));
static constexpr std::double_t X4_const = X1_const * (p_to_3disk/(1-p_disk_to_same_disk));
static constexpr std::double_t X5_const = X1_const * p_to_kdisks + p_disk_to_k_disk * (X2_const + X3_const + X4_const);
static constexpr std::double_t mi1 = 1/Sp;
static constexpr std::double_t mi2 = 1/Sd1;
static constexpr std::double_t mi3 = 1/Sd2;
static constexpr std::double_t mi4 = 1/Sd3;
static constexpr std::double_t mik = 1/Sdk;

std::uint64_t index_max(std::double_t * arr_U, std::uint64_t n)
{
    std::double_t max = arr_U[0];
    std::uint64_t index = 0;
    for(std::uint64_t i = 1; i < n; i++)
        if(max < arr_U[i])
        {
            index = i;
            max = arr_U[i];
        }
    return index;
}
std::double_t U_max(std::double_t * arr_U, std::uint64_t n)
{
    std::double_t max = arr_U[0];
    for(std::uint64_t i = 1; i < n; i++)
        if(max < arr_U[i])
        {
            max = arr_U[i];
        }
    return max;
}

void getSubMatrix(std::double_t** src, std::double_t** dest, std::uint64_t skipRow, std::uint64_t skipCol, std::uint64_t N) 
{
    std::uint64_t i = 0, j = 0;
    for (std::uint64_t row = 0; row < N; row++) 
    {
        for (std::uint64_t col = 0; col < N; col++) 
        {
            if (row != skipRow && col != skipCol) 
            {
                dest[i][j++] = src[row][col];
                if (j == N - 1) 
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

std::double_t getDeterminant(std::double_t** Mat, std::uint64_t N) 
{
    if (N == 1) return Mat[0][0];
    if (N == 2) return (Mat[0][0] * Mat[1][1]) - (Mat[0][1] * Mat[1][0]);

    std::double_t det = 0;
    std::double_t** subMat = new std::double_t*[N - 1];
    for (std::uint64_t i = 0; i < N - 1; i++) 
        subMat[i] = new std::double_t[N - 1];

    std::int8_t sign = 1;
    for (std::uint64_t f = 0; f < N; f++) 
    {
        getSubMatrix(Mat, subMat, 0, f, N);
        det += sign * Mat[0][f] * getDeterminant(subMat, N - 1);
        sign = -sign;
    }

    for (std::uint64_t i = 0; i < N - 1; i++) delete[] subMat[i];
    delete[] subMat;

    return det;
}

std::double_t** getInverseMatrix(std::double_t** Mat, std::uint64_t N) 
{

    std::double_t det = getDeterminant(Mat, N);
    if (std::abs(det) < 1e-12) 
    {
        std::cerr << "Matrix is singular, inverse doesn't exist!" << std::endl;
        return nullptr;
    }

    std::double_t** inverse = new std::double_t*[N];
    
    for (std::uint64_t i = 0; i < N; i++) 
    {
        inverse[i] = new std::double_t[N]();
    }

    if (N == 1) 
    {
        inverse[0][0] = 1.0 / det;
        return inverse;
    }

    std::double_t** subMat = new std::double_t*[N - 1];
    
    for (std::uint64_t i = 0; i < N - 1; i++) 
    {
        subMat[i] = new std::double_t[N - 1];
    }

    for (std::uint64_t i = 0; i < N; i++) 
    {
        for (std::uint64_t j = 0; j < N; j++) 
        {
            getSubMatrix(Mat, subMat, i, j, N);

            std::int8_t sign = ((i + j) % 2 == 0) ? 1 : -1;

            inverse[j][i] = (sign * getDeterminant(subMat, N - 1)) / det;
        }
    }

    for (std::uint64_t i = 0; i < N - 1; i++) delete[] subMat[i];
    delete[] subMat;

    return inverse;
}

std::double_t* multiplyVectorMatrix(std::double_t* vec, std::double_t** mat, std::uint64_t N) 
{
    std::double_t* result = new std::double_t[N]();

    for (std::uint64_t j = 0; j < N; j++) 
    {
        for (std::uint64_t k = 0; k < N; k++) 
        {
            result[j] += vec[k] * mat[k][j];
        }
    }

    return result;
}

void analitical(std::uint64_t r_index , std::uint64_t K_index, std::ofstream& protoci, std::ofstream& grafik, std::ofstream& sve )
{
    std::double_t r = arr_of_r_vals[r_index];
    std::uint64_t K = arr_of_k_vals[K_index];

    std::uint64_t N = K + 4;

    std::double_t** M_mat = new std::double_t* [N]();
    std::double_t** Inverse_M_mat = new std::double_t* [N]();
    std::double_t** P_mat = new std::double_t* [N]();
    std::double_t** I_mat = new std::double_t* [N]();
    std::double_t* A_vec = new std::double_t  [N]();
    
    A_vec[0] = 1;

    for(std::uint64_t i = 0 ; i < N; i++)
    {
        M_mat[i] = new std::double_t [N]();
        Inverse_M_mat[i] = new std::double_t [N]();
        P_mat[i] = new std::double_t [N]();
        I_mat[i] = new std::double_t [N]();
    }

    for(std::uint64_t i = 0 ; i < N; i++)
        for(std::uint64_t j = 0; j < N; j++)
            if(i == j)I_mat[i][j] = 1;

    P_mat[0][0] = p_to_cpu;
    P_mat[0][1] = p_to_1disk;
    P_mat[0][2] = p_to_2disk;
    P_mat[0][3] = p_to_3disk;


    for(std::uint64_t i = 1 ; i < 4; i++)
    {
        P_mat[i][i] = p_disk_to_same_disk;
        P_mat[i][0] = p_disk_to_cpu;
    }

    for(std::uint64_t i = 0 ; i < K; i++)
        P_mat[0][i+4] = p_to_kdisks/K;

    for(std::uint64_t i = 1 ; i < 4; i++)
        for(std::uint64_t j = 0; j < K ; j++)
            P_mat[i][j+4] = p_disk_to_k_disk/K;

    for(std::uint64_t i = 0 ; i < N; i++)
    {
        for(std::uint64_t j = 0 ; j < N; j++)
            M_mat[i][j] = I_mat[i][j] - P_mat[i][j];
    }

    Inverse_M_mat = getInverseMatrix(M_mat, N);

    std::double_t* X = multiplyVectorMatrix(A_vec, Inverse_M_mat, N);

    std::double_t X1_function_of_a = X[0];
    std::double_t X2_function_of_a = X[1];
    std::double_t X3_function_of_a = X[2];
    std::double_t X4_function_of_a = X[3];
    std::double_t X5_function_of_a = X[4] * K;
    std::double_t X5_k_function_of_a = X5_function_of_a / K;

    std::double_t U1_function_of_a = X1_function_of_a / mi1;
    std::double_t U2_function_of_a = X2_function_of_a / mi2;
    std::double_t U3_function_of_a = X3_function_of_a / mi3;
    std::double_t U4_function_of_a = X4_function_of_a / mi4;
    std::double_t U5_k_function_of_a = X5_k_function_of_a / mik;

    std::double_t ro1_function_of_a =   X1_function_of_a    / mi1;
    std::double_t ro2_function_of_a =   X2_function_of_a    / mi2;
    std::double_t ro3_function_of_a =   X3_function_of_a    / mi3;
    std::double_t ro4_function_of_a =   X4_function_of_a    / mi4;
    std::double_t ro5_k_function_of_a = X5_k_function_of_a  / mik;

    std::double_t J1_function_of_a = ro1_function_of_a/(1-ro1_function_of_a);
    std::double_t J2_function_of_a = ro2_function_of_a/(1-ro2_function_of_a);
    std::double_t J3_function_of_a = ro3_function_of_a/(1-ro3_function_of_a);
    std::double_t J4_function_of_a = ro4_function_of_a/(1-ro4_function_of_a);
    std::double_t J5_k_function_of_a = ro5_k_function_of_a/(1-ro5_k_function_of_a);

    std::double_t T1_function_of_a = 1/(mi1 - X1_function_of_a);
    std::double_t T2_function_of_a = 1/(mi2 - X2_function_of_a);
    std::double_t T3_function_of_a = 1/(mi3 - X3_function_of_a);
    std::double_t T4_function_of_a = 1/(mi4 - X4_function_of_a);
    std::double_t T5_k_function_of_a = 1/(mik - X5_k_function_of_a);


    std::double_t R_function_of_a = T1_function_of_a * X1_const +
    T2_function_of_a * X2_const + T3_function_of_a * X3_const +
    T4_function_of_a * X4_const + T5_k_function_of_a * X5_const;
     
    std::double_t U_arr[5] = {U1_function_of_a,U2_function_of_a,U3_function_of_a,U4_function_of_a,U5_k_function_of_a};
    std::double_t Umax = U_max(U_arr, 5);
    std::uint64_t amax_index = index_max(U_arr, N);
    std::double_t amax = (1.0/Umax);
    grafik << r << (r == 0.55 ? "\t\t" :"\t\t\t") << K << "\t" << amax << std::endl;
    amax = amax * r;
    std::double_t U[5] = {U1_function_of_a * amax,U2_function_of_a * amax ,U3_function_of_a * amax,U4_function_of_a * amax,U5_k_function_of_a * amax};
    std::double_t J[5] = {U[0]/(1-U[0]),U[1]/(1-U[1]) ,U[2]/(1-U[2]),U[3]/(1-U[3]),U[4]/(1-U[4])};
    std::double_t T[5] = {1/(mi1 - X[0] * amax), 1/(mi2 - X[1]*amax),1/(mi3 - X[2] * amax), 1/(mi4 - X[3]*amax),1/(mik - X[4] * amax)};

    std::double_t R = T[0] * X1_const +
    T[1] * X2_const + T[2] * X3_const +
    T[3] * X4_const + T[4] * X5_const;


    protoci << "r = " << r << " K = " << K << " Protoci = {";
    for(std::uint64_t i = 0; i < N - 1; i++)
    {
        protoci << " ";
        switch (i) 
        {
            case 0: {protoci << X1_function_of_a;break;}
            case 1: {protoci << X2_function_of_a;break;}
            case 2: {protoci << X3_function_of_a;break;}
            case 3: {protoci << X4_function_of_a;break;}
            case 4: {protoci << X5_k_function_of_a;break;}
            case 5: {protoci << X5_k_function_of_a;break;}
            case 6: {protoci << X5_k_function_of_a;break;}
            case 7: {protoci << X5_k_function_of_a;break;}
        }
        protoci << ",";
    }
    protoci << " " << X5_k_function_of_a << " }" << std::endl;
    for(std::uint64_t i = 0 ; i < 217; i++)protoci<<'-';
    protoci << std::endl;

    

    sve << "r = " << r << " K = " << K << std::endl;
    for(std::uint64_t i = 0 ; i < 217; i++)sve<<'-';
    sve << std::endl;
    sve << "Protoci u fuknciji od a = {";
    for(std::uint64_t i = 0; i < N - 1; i++)
    {
        sve << " ";
        switch (i) 
        {
            case 0: {sve << X1_function_of_a;break;}
            case 1: {sve << X2_function_of_a;break;}
            case 2: {sve << X3_function_of_a;break;}
            case 3: {sve << X4_function_of_a;break;}
            case 4: {sve << X5_k_function_of_a;break;}
            case 5: {sve << X5_k_function_of_a;break;}
            case 6: {sve << X5_k_function_of_a;break;}
            case 7: {sve << X5_k_function_of_a;break;}
        }
        sve << ",";
    }
    sve << " " << X5_k_function_of_a << " } " << std::endl;
    sve << "Protoci za amax = {";
    for(std::uint64_t i = 0; i < N - 1; i++)
    {
        sve << " ";
        switch (i) 
        {
            case 0: {sve << X1_function_of_a * amax;break;}
            case 1: {sve << X2_function_of_a * amax;break;}
            case 2: {sve << X3_function_of_a * amax;break;}
            case 3: {sve << X4_function_of_a * amax;break;}
            case 4: {sve << X5_k_function_of_a * amax;break;}
            case 5: {sve << X5_k_function_of_a * amax;break;}
            case 6: {sve << X5_k_function_of_a * amax;break;}
            case 7: {sve << X5_k_function_of_a * amax;break;}
        }
        sve << ",";
    }
    sve << " " << X5_k_function_of_a * amax << " }" << std::endl;
    
    sve << "Iskoristenje u fuknciji od a = {";
    for(std::uint64_t i = 0; i < N - 1; i++)
    {
        sve << " ";
        switch (i) 
        {
            case 0: {sve << U1_function_of_a;break;}
            case 1: {sve << U2_function_of_a;break;}
            case 2: {sve << U3_function_of_a;break;}
            case 3: {sve << U4_function_of_a;break;}
            case 4: {sve << U5_k_function_of_a;break;}
            case 5: {sve << U5_k_function_of_a;break;}
            case 6: {sve << U5_k_function_of_a;break;}
            case 7: {sve << U5_k_function_of_a;break;}
        }
        sve << ",";
    }
    sve << " " << U5_k_function_of_a << " } " << std::endl;
    sve << "Iskoristenje za amax = {";
    for(std::uint64_t i = 0; i < N - 1; i++)
    {
        sve << " ";
        switch (i) 
        {
            case 0: {sve << U1_function_of_a * amax;break;}
            case 1: {sve << U2_function_of_a * amax;break;}
            case 2: {sve << U3_function_of_a * amax;break;}
            case 3: {sve << U4_function_of_a * amax;break;}
            case 4: {sve << U5_k_function_of_a * amax;break;}
            case 5: {sve << U5_k_function_of_a * amax;break;}
            case 6: {sve << U5_k_function_of_a * amax;break;}
            case 7: {sve << U5_k_function_of_a * amax;break;}
        }
        sve << ",";
    }
    sve << " " << U5_k_function_of_a * amax << " }" << std::endl;
    
    sve << "Prosecan broj poslova u fuknciji od a = {";
    for(std::uint64_t i = 0; i < N - 1; i++)
    {
        sve << " ";
        switch (i) 
        {
            case 0: {sve << J1_function_of_a;break;}
            case 1: {sve << J2_function_of_a;break;}
            case 2: {sve << J3_function_of_a;break;}
            case 3: {sve << J4_function_of_a;break;}
            case 4: {sve << J5_k_function_of_a;break;}
            case 5: {sve << J5_k_function_of_a;break;}
            case 6: {sve << J5_k_function_of_a;break;}
            case 7: {sve << J5_k_function_of_a;break;}
        }
        sve << ",";
    }
    sve << " " << J5_k_function_of_a << " } " << std::endl;
    sve << "Prosecan broj poslova za amax = {";
    for(std::uint64_t i = 0; i < N - 1; i++)
    {
        sve << " ";
        switch (i) 
        {
            case 0: {sve << J[0];break;}
            case 1: {sve << J[1];break;}
            case 2: {sve << J[2];break;}
            case 3: {sve << J[3];break;}
            case 4: {sve << J[4];break;}
            case 5: {sve << J[4];break;}
            case 6: {sve << J[4];break;}
            case 7: {sve << J[4];break;}
        }
        sve << ",";
    }
    sve << " " << J[4] << " }" << std::endl;
    sve << "Vreme odziva u fuknciji od a = " << R_function_of_a << std::endl;
    sve << "Vreme odziva za amax = " << R << std::endl;
    sve << (amax_index != 4 ? "Kriticni resurs je " : "Kriticni resursi su ") << (amax_index == 0 ? "procesor": (amax_index != 4 ? "sistemski disk" : "korisnici diskovi"));
    if(amax_index != 0 && amax_index != 4) sve << amax_index;
    sve << std::endl;
    
    sve << "Srednje vreme odziva u funkciji od a = " << R_function_of_a << std::endl;
    sve << "Srednje vreme odziva za amax = " << R << std::endl;
    
    for(std::uint64_t i = 0 ; i < 217; i++)sve<<'-';
    sve << std::endl;

    for(std::uint64_t i = 0; i < N;i++)
    {
        delete [] M_mat[i];
        delete [] Inverse_M_mat[i];
        delete [] P_mat[i];
        delete [] I_mat[i];
    }

    delete [] M_mat;
    delete [] Inverse_M_mat;
    delete [] P_mat;
    delete [] I_mat;
    delete [] X;
    delete [] A_vec;
}

std::double_t random_prob()
{
    return (static_cast<std::double_t>(rand()))/(static_cast<std::double_t>(RAND_MAX)) ;
}

static constexpr std::double_t p_to_disk1_sim =  p_to_1disk + p_to_cpu;
static constexpr std::double_t p_to_disk2_sim =  p_to_2disk + p_to_disk1_sim;
static constexpr std::double_t p_to_disk3_sim =  p_to_3disk + p_to_disk2_sim;
static constexpr std::double_t p_to_disksk_sim =  p_to_kdisks + p_to_disk2_sim;
static constexpr std::double_t p_disk_to_same_disk_sim = p_disk_to_cpu + p_disk_to_same_disk;
static constexpr std::double_t to_get_k = 1 - p_disk_to_same_disk_sim;
std::uint64_t prob_to_k(std::double_t p, std::uint64_t K)
{
    std::double_t temp = p - p_disk_to_same_disk_sim;
    switch (K) 
    {
        case 2:
        {
            if(temp <= to_get_k/2)return 0;
            else return 1;
        }
        case 3:
        {
            if(temp <= to_get_k/3)return 0;
            else if(temp <= to_get_k*2/3)return 1;
            else return 2;
        }
        case 4:
        {
            if(temp <= to_get_k/4)return 0;
            else if(temp <= to_get_k*2/4)return 1;
            else if(temp <= to_get_k*3/4)return 2;
            else return 3;
        }
        case 5:
        {
            if(temp <= to_get_k/5)return 0;
            else if(temp <= to_get_k*2/5)return 1;
            else if(temp <= to_get_k*3/5)return 2;
            else if(temp <= to_get_k*4/5)return 3;
            else return 4;
        }
        default:{
            exit(0);
        }
    }
}

struct job
{
    std::uint64_t next = 0;
    std::double_t prev_prob = 0;
    std::uint64_t index;
    std::double_t arrival_time;
    bool end = false;
    job(std::double_t curr_time,std::uint64_t index): index(index),
      arrival_time(curr_time){}
    void to_next(std::uint64_t K)
    {
        switch (next) 
        {

            case 0:
            {
                std::double_t p_next = random_prob();
                prev_prob = p_next; 
                if(p_next <= p_to_cpu)next = 0;
                else if(p_to_cpu < p_next && p_next <= p_to_disk1_sim)next = 1;
                else if(p_to_cpu < p_next && p_next <= p_to_disk2_sim)next = 2;
                else if(p_to_cpu < p_next && p_next <= p_to_disk3_sim)next = 3;
                else next = 4 + prob_to_k(prev_prob,K);
                break;
            }
            case 1:
            {
                std::double_t p_next = random_prob();
                prev_prob = p_next;
                if(p_next <= p_disk_to_cpu)next = 0;
                else if(p_disk_to_cpu < p_next && p_next <= p_disk_to_same_disk_sim)next = 1;
                else next = 4 + prob_to_k(prev_prob,K);;
                break;
            }
            case 2:
            {
                std::double_t p_next = random_prob();
                prev_prob = p_next;
                if(p_next <= p_disk_to_cpu)next = 0;
                else if(p_disk_to_cpu < p_next && p_next <= p_disk_to_same_disk_sim)next = 2;
                else next = 4 + prob_to_k(prev_prob,K);;
                break;
            }
            case 3:
            {
                std::double_t p_next = random_prob();
                prev_prob = p_next;
                if(p_next <= p_disk_to_cpu)next = 0;
                else if(p_disk_to_cpu < p_next && p_next <= p_disk_to_same_disk_sim)next = 3;
                else next = 4 + prob_to_k(prev_prob,K);
                break;
            }
            default:
            {
                end = true;
            }
        }
    }
};

struct queue_data
{
    std::double_t time_worked = 0;
    std::double_t wait_time = 0;
    std::double_t time_to_process;
    std::uint64_t job_index = 0;
    bool end = false;
    queue_data(std::double_t time_to_process, std::uint64_t index):time_to_process(time_to_process),job_index(index){};
    std::double_t process(std::double_t time)
    {
        std::double_t temp = time_worked;
        time_worked = time_to_process > (time + time_worked) ? time + time_worked : 0;
        end  = time_worked == 0 ? true : false;
        return  time_worked == 0 ? time - (time_to_process - temp) : 0;
    }
    void wait(std::double_t time)
    {
        wait_time += time;
    }
};

struct comp
{
    std::double_t time_on = 0;
    std::uint64_t finished = 0;
    std::uint64_t waiting = 0;
    std::double_t total_time = 0;
    bool on = false;
};

std::double_t exponential(std::double_t lambda)
{
    return -log(1.0 - ( ( (double_t)(uint64_t)rand()) + 1 ) / ((uint64_t)RAND_MAX + 2) ) / lambda;
}

std::double_t ret_index_to_time(std::uint64_t index)
{
    switch (index) 
    {
    case 0: return Sp;
    case 1: return Sd1;
    case 2: return Sd2;
    case 3: return Sd3;
    default: return Sdk;
    }
}
std::double_t total_time = 0.0;
std::uint64_t completed = 0;
comp* simulation(std::uint64_t r_index , std::uint64_t K_index, std::uint64_t vreme = 30)
{
    std::uint64_t K     = arr_of_k_vals[K_index];
    std::double_t r     = arr_of_r_vals[r_index];
    std::double_t timeout   = vreme * 60;
    std::uint64_t N     = K + 4;
    std::double_t U_a[5] = {X1_const/mi1,X2_const/mi2,X3_const/mi3,X4_const/mi4,(X5_const/K)/mik};
    std::double_t a = r*(1.0 / U_max(U_a, N));

    comp* ret= new comp [9](); 
    std::double_t curr_time = 0.0;
    std::double_t prev_time = 0.0;
    std::double_t next_arrival = exponential(a);
    std::vector<job> jobs;
    std::vector<queue_data> queues[9] = 
    {std::vector<queue_data>(),std::vector<queue_data>(),std::vector<queue_data>(),
    std::vector<queue_data>(),std::vector<queue_data>(),std::vector<queue_data>(),
    std::vector<queue_data>(),std::vector<queue_data>(),std::vector<queue_data>()};
    std::double_t next_time = Sdk;
    while (curr_time <= timeout) 
    {

        std::double_t passed_time = curr_time - prev_time;

        for(std::uint64_t i = 0 ; i < 9 ; i++)
        {
            if(ret[i].on == false)continue;
            std::double_t loop_time = passed_time;
            while(loop_time > 0)
            {
                std::double_t rem_time = queues[i][0].time_to_process - queues[i][0].time_worked;
                std::double_t temp = queues[i][0].process( loop_time);
                if(queues[i][0].end)
                {
                    ret[i].waiting -= 1;
                    ret[i].finished += 1;
                    if(ret[i].waiting == 0)ret[i].on = false;
                    std::uint64_t ind = queues[i][0].job_index;
                    jobs[ind].to_next(K);
                    if(jobs[ind].end == false)
                    {
                        std::uint64_t ret_index = jobs[ind].next;;
                        ret[ret_index].waiting += 1;
                        if(ret[ret_index].waiting == 1)
                        {
                            ret[ret_index].on = true;
                        }
                        queues[ret_index].push_back(queue_data(ret_index_to_time(ret_index),ind));
                    }
                    else {
                    {
                        total_time += (curr_time - jobs[ind].arrival_time);
                        completed += 1;
                    }
                    }
                    queues[i].erase(queues[i].begin());
                    if(queues[i].size() == 0)
                    {
                            ret[i].on = false;
                            loop_time -= rem_time;
                            break;
                    }for(auto& a : queues[i])
                    {
                        a.wait(loop_time - temp);
                    }
                }
                else
                {
                    for(size_t j=1;j<queues[i].size();j++)
                        queues[i][j].wait(loop_time - temp);
                }

                loop_time = temp;
            }
            ret[i].time_on += passed_time - loop_time;
        }

        while(curr_time >= next_arrival)
        {
                jobs.push_back(job(curr_time,jobs.size()));
                ret[0].waiting += 1;
                queues[0].push_back(queue_data(Sp,jobs.back().index));
                if(ret[0].waiting == 1)
                {
                    ret[0].on = true;
                }
                next_arrival += exponential(a);
        }
        if(next_arrival - curr_time < next_time)
        {
            next_time = next_arrival - curr_time;
        }
        
        for(std::uint64_t i = 0 ; i < 9; i++)
        {
            if(ret[i].on == false)continue;
            if(queues[i][0].time_to_process - queues[i][0].time_worked < next_time)
                next_time = queues[i][0].time_to_process - queues[i][0].time_worked;
        }
        if(next_time <= 1e-12)
            next_time = 1e-12;
        prev_time = curr_time;
        curr_time += next_time;
        next_time = Sdk;
    }
    return ret;
}   

std::double_t U_global[9] = {0,0,0,0,0,0,0,0,0};
std::double_t J_global[9] = {0,0,0,0,0,0,0,0,0};
std::double_t X_global[9] = {0,0,0,0,0,0,0,0,0};
std::double_t T_global[9] = {0,0,0,0,0,0,0,0,0};
std::double_t R_global = 0.0;
void simul_file_write(std::double_t r, std::uint64_t K, comp* ret, std::double_t time, std::ofstream& sim1)
{
        std::uint64_t N = K + 4;

    std::double_t U[9] = {0};
    std::double_t J[9] = {0};
    std::double_t X[9] = {0};
    std::double_t T[9] = {0};
    std::double_t R = 0;
    
    for(std::uint64_t i = 0 ; i < N ; i++)
    {
        U[i] = ret[i].time_on / time;
        
        X[i] = static_cast<std::double_t>(ret[i].finished) / time;
        
        if (ret[i].finished > 0) 
        {
            T[i] = ret[i].time_on / ret[i].finished; 
            J[i] = X[i] * T[i]; 
        } 
        else 
        {
            T[i] = 0.0;
            J[i] = 0.0;
        }
    }

    R = completed > 0 ? (total_time / completed) : 0.0;
    std::uint64_t indexmax = index_max(U, N);

    sim1 << "r = " << r << " K = " << K << std::endl;
    for(std::uint64_t i = 0 ; i < 217; i++)sim1<<'-';
    sim1 << std::endl;
    sim1 << "X = {";
    for(std::uint64_t i = 0; i < N-1 ; i++)
    {

        sim1 << " " << X[i] << ",";
    }
    sim1 << " " << X[N-1] << "}" << std::endl;
    
    sim1 << "U = {";
    for(std::uint64_t i = 0; i < N-1 ; i++)
    {

        sim1 << " " << U[i] << ",";
    }
        sim1 << " " << U[N-1] << "}" << std::endl;

    sim1 << "J = {";
    for(std::uint64_t i = 0; i < N-1 ; i++)
    {

        sim1 << " " << J[i] << ",";
    }
        sim1 << " " << J[N-1] << "}" << std::endl;

    sim1 << "T = {";
    for(std::uint64_t i = 0; i < N-1 ; i++)
    {

        sim1 << " " << T[i] << ",";
    }
        sim1 << " " << T[N-1] << "}" << std::endl;
    sim1 << "R = " << R << std::endl;


    sim1 << "Kriticni resurs je "; 
    if(indexmax == 0) sim1 << "procesor";
    else if (indexmax > 0 && indexmax < 4) sim1 << "sistemski disk " << indexmax;
    else sim1 << "korisicki disk " << indexmax - 4;
    sim1 << std::endl;

    for(std::uint64_t i = 0 ; i < 217; i++)sim1<<'-';
    sim1 << std::endl;


    for(std::uint64_t i = 0; i < N; i++)
    {
        U_global[i] += U[i];
        X_global[i] += X[i];
        J_global[i] += J[i];
        T_global[i] += T[i];
    }
    R_global += R;

}
void simul_calc(std::double_t r, std::uint64_t K, comp* ret, std::double_t time)
{
    std::uint64_t N = K + 4;

    std::double_t U[9] = {0};
    std::double_t J[9] = {0};
    std::double_t X[9] = {0};
    std::double_t T[9] = {0};
    std::double_t R = 0;
    
    for(std::uint64_t i = 0 ; i < N ; i++)
    {
        U[i] = ret[i].time_on / time;
        
        X[i] = static_cast<std::double_t>(ret[i].finished) / time;
        
        if (ret[i].finished > 0) 
        {
            T[i] = ret[i].time_on / ret[i].finished; 
            J[i] = X[i] * T[i]; 
        } 
        else 
        {
            T[i] = 0.0;
            J[i] = 0.0;
        }
    }

    R = completed > 0 ? (total_time / completed) : 0.0;

    for(std::uint64_t i = 0; i < N; i++)
    {
        U_global[i] += U[i];
        X_global[i] += X[i];
        J_global[i] += J[i];
        T_global[i] += T[i];
    }
    R_global += R;
}
int main()
{
    std::ofstream protoci("protoci_analiticki.txt");
    std::ofstream grafik("amax.txt");
    std::ofstream sve("rezultati_analiticki.txt");
    std::ofstream sim1("rezultati_simulacija.txt");
    std::ofstream sim2("rezultati_simulacija_usrednjeno.txt");
    grafik << "r\t\t\tK\tamax" << std::endl;
    analitical(0, 0, protoci, grafik, sve);
    analitical(0, 1, protoci, grafik, sve);
    analitical(0, 2, protoci, grafik, sve);
    analitical(0, 3, protoci, grafik, sve);
    analitical(1, 0, protoci, grafik, sve);
    analitical(1, 1, protoci, grafik, sve);
    analitical(1, 2, protoci, grafik, sve);
    analitical(1, 3, protoci, grafik, sve);
    analitical(2, 0, protoci, grafik, sve);
    analitical(2, 1, protoci, grafik, sve);
    analitical(2, 2, protoci, grafik, sve);
    analitical(2, 3, protoci, grafik, sve);
    analitical(3, 0, protoci, grafik, sve);
    analitical(3, 1, protoci, grafik, sve);
    analitical(3, 2, protoci, grafik, sve);
    analitical(3, 3, protoci, grafik, sve);
    
    comp* ret = new comp[9];
    std::double_t time = 0;
    std::double_t simul_time = 30;

    for(std::uint64_t i = 0 ; i < 4; i++)
        for(std::uint64_t j = 0 ; j < 4; j++)
            {
                total_time = 0;
                completed = 0;
                comp* temp = simulation(i, j,simul_time);
                simul_file_write(arr_of_r_vals[i], arr_of_k_vals[j], temp, simul_time * 60, sim1);
                delete [] temp;
            }

    for(std::uint64_t i = 0 ; i < 4; i++)
        for(std::uint64_t j = 0 ; j < 4; j++)
            {
                for(std::uint64_t k = 0; k < 9 ; k++)
                    U_global[k] = J_global[k] = T_global[k] = X_global[k] = 0.0;
                R_global = 0.0;
                for(std::uint64_t k = 0 ; k < 100; k++)
                {
                    total_time = 0;
                    completed = 0;
                    comp* temp = simulation(i, j,simul_time);
                    simul_calc(arr_of_r_vals[i], arr_of_k_vals[j], temp, simul_time * 60);
                    delete [] temp;
                }
                for(std::uint64_t k = 0; k < 9 ; k++)
                {
                    U_global[k] /= 100;
                    J_global[k] /= 100;
                    T_global[k] /= 100; 
                    X_global[k] /= 100;
                }
                R_global /= 100;
                std::uint64_t indexmax = index_max(U_global, 4 + arr_of_k_vals[j]);
                sim2 << "r = " << arr_of_r_vals[i] << " K = " << arr_of_k_vals[j] << std::endl;
                    for(std::uint64_t i = 0 ; i < 217; i++)sim2<<'-';
                    sim2 << std::endl;
                std::uint64_t N = arr_of_k_vals[j] + 4;
                sim2 << "X = {";
                for(std::uint64_t i = 0; i < N-1 ; i++)
                {

                    sim2 << " " << X_global[i] << ",";
                }
                sim2 << " " << X_global[N-1] << "}" << std::endl;
                
                sim2 << "U = {";
                for(std::uint64_t i = 0; i < N-1 ; i++)
                {

                    sim2 << " " << U_global[i] << ",";
                }
                    sim2 << " " << U_global[N-1] << "}" << std::endl;

                sim2 << "J = {";
                for(std::uint64_t i = 0; i < N-1 ; i++)
                {

                    sim2 << " " << J_global[i] << ",";
                }
                    sim2 << " " << J_global[N-1] << "}" << std::endl;

                sim2 << "T = {";
                for(std::uint64_t i = 0; i < N-1 ; i++)
                {

                    sim2 << " " << T_global[i] << ",";
                }
                sim2 << " " << T_global[N-1] << "}" << std::endl;
                
                sim2 << "R = " << R_global << std::endl;
                
                sim2 << "Kriticni resurs je "; 
                if(indexmax == 0) sim2 << "procesor";
                else if (indexmax > 0 && indexmax < 4) sim2 << "sistemski disk " << indexmax;
                else sim2 << "korisicki disk " << indexmax - 4;
                sim2 << std::endl;
                
                for(std::uint64_t i = 0 ; i < 217; i++)sim2<<'-';
                    sim2 << std::endl;
            }
    
    sim1.close();
    sim2.close();
    protoci.close();
    grafik.close();
    sve.close();
    return 0;
}
