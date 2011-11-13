/* C source code
 * File: "/home/kassick/Work/olam/trace2paje/inputs/filenames/filenames.c"
 * Created: "Sáb, 12 Nov 2011 22:55:46 -0200 (kassick)"
 * Updated: "Dom, 13 Nov 2011 03:46:21 -0200 (kassick)"
 * $Id$
 * Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
 */
/*
 * ===========================================================================
 *
 *       Filename:  teste.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  12-11-2011 22:55:46 BRST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include "rst_functions.h"
#include "evt_defs.h"


#define TEST_SIZE 100000
#define FILE_SIZE 10000
int main(int argc, char *argv[])
{
  long * avector;
  int i,cur_vect_pos, cur_file;
  int * fdlist;
  char fname[100];


  rst_init(0,0);

  // fill in some data
  avector = malloc(sizeof(long) * TEST_SIZE);
  fdlist = malloc(sizeof(int) * (TEST_SIZE/FILE_SIZE) + 1);

  for (i = 0; i < TEST_SIZE; i++)
    avector[i] = i*i; // don't care if overflow

  cur_vect_pos = 0;
  cur_file = 0;
  while (cur_vect_pos < TEST_SIZE)
  {
    sprintf(fname,"file_%d.out",cur_file);
    rst_event_s(OPEN_IN,fname);
    fdlist[cur_file] = open(fname, O_CREAT | O_WRONLY, S_IRWXU);
    if (fdlist[cur_file] == -1)
    {
      fprintf(stderr,"Can not open a file %s : %s\n",fname, strerror(errno));
      exit(1);
    }
    rst_event_is(OPEN_OUT,fdlist[cur_file], fname);

    rst_event_i(WRITE_IN,fdlist[cur_file]);
    write(fdlist[cur_file], 
        avector + cur_vect_pos, 
        FILE_SIZE*sizeof(long));
    rst_event_i(WRITE_OUT,fdlist[cur_file]);


    cur_vect_pos += FILE_SIZE;
    cur_file++;
  }

  cur_file --;

  while (cur_file >= 0)
  {
    rst_event_i(CLOSE_IN,fdlist[cur_file]);
    close(fdlist[cur_file]);
    rst_event_i(CLOSE_OUT,fdlist[cur_file]);
    cur_file--;
  }

  rst_finalize();

  return 0;
}
