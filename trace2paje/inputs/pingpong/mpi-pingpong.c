/* C source code
 * File: "/home/kassick/Work/olam/trace2paje/inputs/pingpong/mpi-pingpong.c"
 * Created: "Ter, 04 Out 2011 17:13:06 -0300 (kassick)"
 * Updated: "Qua, 05 Out 2011 14:29:02 -0300 (kassick)"
 * $Id$
 * Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
 */
/*
 * ===========================================================================
 *
 *       Filename:  mpi-pingpong.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04-10-2011 17:13:06 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include <mpi.h>
#include <stdio.h>
#include "evt_defs.h"
#include <rastro.h>
#include "rst_functions.h"

#define MAX_HOSTNAME 200
#define SRC 0
#define DST 1

int main(int argc, char ** argv)
{
  char hostname[MAX_HOSTNAME];
  int hostname_len = MAX_HOSTNAME;
  int rank;

  int tag = 0;
  int data = 1;
  
  MPI_Status status;




  MPI_Init(&argc, &argv);

  MPI_Get_processor_name(hostname,&hostname_len);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  rst_init(rank,10); // 2nd is ignore

  rst_event_s(INIT_N, hostname);

  for(tag = 0; tag < 100; tag++)
  {
    if (rank == SRC)
    {
      rst_event_ii(PING_IN,tag,DST);

      MPI_Send(&data,1,MPI_INT, DST, tag, MPI_COMM_WORLD);
      printf("Rank %d sent %d, tag %d\n",rank,data,tag);
      tag++;
      data+=10;
    } else {
      MPI_Recv(&data,1,MPI_INT, SRC, tag, MPI_COMM_WORLD,&status);
      
      rst_event_ii(PING_OUT,tag,DST);

      printf("Rank %d got %d, tag %d\n",rank,data,tag);
      tag++;
      data+=10;
    }
  }


  MPI_Finalize();
  rst_event_s(FINALIZE_N, hostname);

  rst_finalize();


  return 0;
}
