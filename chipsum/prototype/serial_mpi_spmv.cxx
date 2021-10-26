/*
 * @Author: Li Kunyun
 * @Date: 2021-09-24 14:59:41
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-09-24 15:46:40
 * @Description: 
 */
#include <mpi.h>
#include <vector>



struct MPIPack
{

    std::vector<int> elements_to_send;
    std::vector<int>            neighbors;
    std::vector<int>  recv_length;
    std::vector<int>  send_length;
    std::vector<double>        send_buffer;
    std::vector<MPI_Request>   request;
};





void ExchangeExternels(MPIPack& A,
                  std::vector<double>& x,int local_nrow)
{


  int numprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  if (numprocs < 2) return;

  
  int num_neighbors = A.neighbors.size();
  const std::vector<int>& recv_length = A.recv_length;
  const std::vector<int>& send_length = A.send_length;
  const std::vector<int>& neighbors = A.neighbors;
  const std::vector<int>& elements_to_send = A.elements_to_send;

  std::vector<double>& send_buffer = A.send_buffer;

  int MPI_MY_TAG = 99;

  std::vector<MPI_Request>& request = A.request;

  
  double* x_external = &(x[local_nrow]);

  MPI_Datatype mpi_dtype = MPI_DOUBLE;

  for(int i=0; i<num_neighbors; ++i) {
    int n_recv = recv_length[i];
    MPI_Irecv(x_external, n_recv, mpi_dtype, neighbors[i], MPI_MY_TAG,
              MPI_COMM_WORLD, &request[i]);
    x_external += n_recv;
  }

  size_t total_to_be_sent = elements_to_send.size();

  for(size_t i=0; i<total_to_be_sent; ++i) {

    send_buffer[i] = x[elements_to_send[i]];
  }

  double* s_buffer = &send_buffer[0];

  for(int i=0; i<num_neighbors; ++i) {
    int n_send = send_length[i];
    MPI_Send(s_buffer, n_send, mpi_dtype, neighbors[i], MPI_MY_TAG,
             MPI_COMM_WORLD);
    s_buffer += n_send;
  }



  MPI_Status status;
  for(int i=0; i<num_neighbors; ++i) {
    if (MPI_Wait(&request[i], &status) != MPI_SUCCESS) {
      std::cerr << "MPI_Wait error\n"<<std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

}

int main(int argc,char* argv[])
{
  /**
   *    
   *    o————o————o——|——o————o 
   *    1    2    3     4    5
   *  
   *  | 1  1  0  0  0 |   | 1 |
   *  | 1  1  1  0  0 |   | 2 |
   *  | 0  1 -1 *1  0 |   | 3 |
   * ————————————————————————————
   *  | 0  0 *1 -1  1 |   | 4 |
   *  | 0  0  0  1  1 |   | 5 |
   * 
   */
  
  MPI_Init(&argc,&argv);

  int id,size;

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&id);

  std::vector<double> x;
 

  MPIPack A;

  int local_nrow ;

  if(id == 0){
    x.resize(4);
    for(int i=0;i<3;++i)x[i] = double(i+1);
    A.elements_to_send.push_back(2);
    A.neighbors.push_back(1);
    A.send_length.push_back(1);
    A.recv_length.push_back(1);
    A.send_buffer.resize(1);
    A.request.resize(1);
    local_nrow = 3;
  }
  else if(id == 1){
    x.resize(3);
    for(int i=0;i<2;++i)x[i] = double(i+4);
    A.elements_to_send.push_back(0);
    A.neighbors.push_back(0);
    A.send_length.push_back(1);
    A.recv_length.push_back(1);
    A.send_buffer.resize(1);
    A.request.resize(1);
    local_nrow = 2;
  }
  
  ExchangeExternels(A,x,local_nrow);


  for(int i=0;i<1000000*id;++i);
  std::cout<<"Proc "<<id<<": ";
  for(int i=0;i<x.size();++i){
    std::cout<<x[i]<<" ";
  }
  std::cout<<std::endl;


  MPI_Finalize();
}