/* C source code
 * File: "/home/kassick/Work/olam/trace2paje/inputs/pingpong-2/mpi-pingpong.c"
 * Created: "Ter, 04 Out 2011 17:13:06 -0300 (kassick)"
 * Updated: "Seg, 10 Out 2011 19:09:13 -0300 (kassick)"
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
#define NELEM 100000000

int main(int argc, char ** argv)
{
  char hostname[MAX_HOSTNAME];
  int hostname_len = MAX_HOSTNAME;
  int rank;
  int * buf;

  int tag = 0;
  
  MPI_Status status;




  MPI_Init(&argc, &argv);

  MPI_Get_processor_name(hostname,&hostname_len);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  rst_init(rank,10); // 2nd is ignore

  rst_event_s(INIT_N, hostname);

  buf = malloc(sizeof(int)*NELEM);
  if (!buf)
  {
    printf("Error! buf is null\n");
    exit(1);
  }

  for(tag = 0; tag < 10; tag++)
  {
    if ((rank % 2) == 0)
    {
      buf[0] = rank;
      buf[NELEM-1] = rank;
      // even rank sends
      rst_event_iii(PING_IN,tag,rank,rank+1);

      MPI_Send(buf,NELEM,MPI_INT, rank+1, tag, MPI_COMM_WORLD);

      rst_event(PING_S_OUT);

      sleep(1);

      printf("Rank %d sent %d to %d, tag %d\n",rank,buf[0],rank+1,tag);
      tag++;
      buf[0]+=10;
    } else {
      rst_event_ii(PING_R_IN,tag,rank-1);
      MPI_Recv(buf,NELEM,MPI_INT, rank-1, tag, MPI_COMM_WORLD,&status);
      
      rst_event_iii(PING_OUT,tag,rank-1,rank);

      printf("Rank %d got %d from %d, tag %d\n",rank,buf[0],rank-1,tag);
      tag++;
    }
  }


  MPI_Finalize();
  rst_event_s(FINALIZE_N, hostname);

  rst_finalize();


  return 0;
}
