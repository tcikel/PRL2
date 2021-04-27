#include <vector>
#include <mpi.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <time.h>


#define MAT1_ROW 0
#define MAT1_COL 1
#define MAT2_ROW 2
#define MAT2_COL 3
#define COLUMN 4
#define ROW 5
#define FINAL 6
#define TEST true
using namespace std;



/*
 * This function is taken from http://www.guyrutenberg.com/2007/09/22/profiling-code-using-clock_gettime/
 * 
*/
timespec diff(timespec start, timespec end)
{
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}


/**
 * @param y y dimension
 * @param x x dimensio
 * @param column number of column in matrix
 * @return 1D index
 */
int getMatrixIndex(int y, int x, int column)
{
    return (y * column) + x;
}

/**
 * Wait for data from other procesors
 * @param recv_from wait on data form this procesor id
 * @param store_data  store received data
 * @param tag either ROW or COLUMN
 * @param received_count
 */
void recv_data(int recv_from,int *store_data, int tag,int *received_count)
{ 
    MPI_Status status;
    MPI_Recv(store_data, 1, MPI_INT, recv_from, tag, MPI_COMM_WORLD, &status);
    (*received_count)++;
    return;
}

/**
 * Main algorithm for all procesors in matrix
 * Receive data, make the calculation and send the data to another procesor
 * @param my_id id of procesor
 */
void mesh_algorithm(int my_id)
{
    int mat1_row, mat1_col,mat2_row,mat2_col = 0;
    int received = 0;
    int c = 0;
    int recv_x, recv_y = 0;
    int data_x, data_y = 0;
    int x, y;
    MPI_Status status;



    //Get dimensions from first procesor
    MPI_Recv(&mat1_row, 1, MPI_INT, 0, MAT1_ROW, MPI_COMM_WORLD, &status);
    MPI_Recv(&mat1_col, 1, MPI_INT, 0, MAT1_COL, MPI_COMM_WORLD, &status);
    MPI_Recv(&mat2_row, 1, MPI_INT, 0, MAT2_ROW, MPI_COMM_WORLD, &status);
    MPI_Recv(&mat2_col, 1, MPI_INT, 0, MAT2_COL, MPI_COMM_WORLD, &status);

    //Calculate my x,y index and index form who i will be receiving
    x = my_id % mat2_col;
    y = my_id / mat2_col;
    recv_x = (x < 1) ? 0 : getMatrixIndex(y, x-1, mat2_col);
    recv_y = (y < 1) ? 0 : getMatrixIndex(y-1, x, mat2_col);

    while (received < (mat1_col + mat2_row))
    {
        //Wait for data
        recv_data(recv_x,&data_x,ROW,&received);
        recv_data(recv_y,&data_y,COLUMN,&received);
        c += data_x * data_y;
        
        //If id is not at border of mesh send to another processor 
        if ( x!=(mat2_col-1) )
        {
            MPI_Send(&data_x, 1, MPI_INT, getMatrixIndex(y, x+1, mat2_col), ROW, MPI_COMM_WORLD);
        } 
        if ( y!=(mat1_row-1) )
        { 
            MPI_Send(&data_y, 1, MPI_INT,  getMatrixIndex(y+1, x, mat2_col), COLUMN, MPI_COMM_WORLD);
        }
    }
    //Send final data to first procesor
    MPI_Send(&c, 1, MPI_INT, 0, FINAL, MPI_COMM_WORLD);
    return;
}



int getMatrix(int *size, vector<int> *matrix, string filename)
{
    string line ;
    int tmp_num;
    
    ifstream matrixFile(filename.c_str());
    getline(matrixFile, line );
    *size = stoi(line );
    
    while (getline(matrixFile, line ))
    {
        istringstream lineStream(line );
        while (lineStream >> tmp_num)
        {
        matrix->push_back(tmp_num);
        }
    }
    return 0;
}


int main(int argc, char *argv[])
{

    timespec start_time;
    timespec end_time;
    int proces_num;
    int my_id = 0;
    vector<int> output;
    int mat1_row;
    int mat1_col;
    int mat2_row;
    int mat2_col;

    //Init OpenMPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proces_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Status status;

    if (my_id == 0)
    {

        vector<int> matrix1;
        vector<int> matrix2;
        getMatrix(&mat1_row, &matrix1, "mat1");
        getMatrix(&mat2_col, &matrix2, "mat2");

        mat1_col = matrix1.size() / mat1_row;
        mat2_row = matrix2.size() / mat2_col;

        //Send matrix dimension to everybody
        for (int i = 0; i < proces_num; i++)
        {
            MPI_Send(&mat1_col, 1, MPI_INT, i, MAT1_COL, MPI_COMM_WORLD);
            MPI_Send(&mat1_row, 1, MPI_INT, i, MAT1_ROW, MPI_COMM_WORLD);
            MPI_Send(&mat2_col, 1, MPI_INT, i, MAT2_COL, MPI_COMM_WORLD);
            MPI_Send(&mat2_row, 1, MPI_INT, i, MAT2_ROW, MPI_COMM_WORLD);
        }

        //Senda data from first matrix to specific procesors
        for (int i = 0; i < mat1_row; i++)
        {
            for (int j = 0; j < mat1_col; j++)
            {
                int tmp_index = getMatrixIndex(i, j, mat1_col);
                int send_id = getMatrixIndex(i, 0, mat2_col);
                MPI_Send(&matrix1[tmp_index], 1, MPI_INT, send_id, ROW, MPI_COMM_WORLD);
            }
        }
        //Send data from second matrix to specific procesors
        for (int i = 0; i < mat2_col; i++)
        {
            for (int j = 0; j < mat2_row; j++)
            {
                int tmp_index = getMatrixIndex(j, i, mat2_col);
                int send_id = getMatrixIndex(0, i, mat2_col);
                MPI_Send(&matrix2[tmp_index], 1, MPI_INT, send_id, COLUMN, MPI_COMM_WORLD);
            }
        }
    }

    if(TEST)
    { 
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
    } 

    //Core of mesh algorithm 
    mesh_algorithm(my_id);

    //First process receive the data from everobdy else and prints it
    if (my_id == 0)
    {
        int tmp;
        for (int i = 0; i < proces_num; i++)
        {
            MPI_Recv(&tmp, 1, MPI_INT, i, FINAL, MPI_COMM_WORLD, &status);
            output.push_back(tmp);
        }
        tmp=0;
        cout << mat1_row << ":" << mat2_col <<endl; 
        for(int i= 0; i<output.size(); i++)
        { 
            tmp++;
            if (tmp==mat2_col)
            { 
                cout<<output[i]<<endl;
                tmp=0;
            }
            else
            {
            cout<<output[i]<<" ";
            }
        }
    }
    //Wait for everbody and calculate time
    MPI_Barrier(MPI_COMM_WORLD);    
    if(TEST && my_id==0){
	    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_time);
        timespec difference = diff(start_time,end_time);
        cout << difference.tv_sec << ":" << difference.tv_nsec <<endl;
    }

    MPI_Finalize();
    return 0;
}
