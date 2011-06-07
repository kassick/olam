/* C source code
 * File: "/home/kassick/Work/olam/olamfs-trace-converter/src/gen_fortran_names.c"
 * Created: "Ter, 07 Jun 2011 18:48:10 -0300 (kassick)"
 * Updated: "Ter, 07 Jun 2011 19:28:00 -0300 (kassick)"
 * $Id$
 * Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
 */
/*
 * ===========================================================================
 *
 *       Filename:  gen_fortran_names.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07-06-2011 18:48:10 BRT
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
#include <getopt.h>
#include <aky.h>
#include <olam_events.h>



#define HEADER "Module rastro_evts\n"
#define FOOTER "End Module rastro_evts\n"

int main(int argc, char **argv)
{
  char *ofname = "/dev/stdout";
  int opt;
  FILE *of;
  int i;

  while ((opt = getopt(argc,argv,"o:")) != -1) {
    switch (opt) {
      case 'o':
        ofname = optarg;
        break;

      default:
        printf("Usage: %s -o filename\n",argv[0]);
        exit(1);
    }
  }

  name_init();
  of = fopen(ofname,"w+");
  if (!of) {
    printf("Can not open output file %s\n",ofname);
    exit(1);
  }

  fprintf(of,HEADER);
  for (i = 0; olam_evt_names[i].name != NULL; i++)
  {
    if (olam_evt_names[i].type == EVT) {
      fprintf(of,"integer :: OLAM_EVT_%s  =  %d\n",olam_evt_names[i].short_name,olam_evt_names[i].id);
    } else if (olam_evt_names[i].type == IN) {
      fprintf(of,"integer :: OLAM_EVT_%s_IN  =  %d\n",olam_evt_names[i].short_name,olam_evt_names[i].id);
    } else {
      fprintf(of,"integer :: OLAM_EVT_%s_OUT  =  %d\n",olam_evt_names[i].short_name,olam_evt_names[i].id);
    }
  }

  fprintf(of,FOOTER);

  fclose(of);

}


