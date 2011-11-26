/* C source code
 * File: "/home/kassick/Work/olam/orangefs-rastropatch/src/libRastro/uniq_ids_q.c"
 * Created: "Dom, 13 Nov 2011 20:56:10 -0200 (kassick)"
 * Updated: "Sáb, 26 Nov 2011 00:36:43 -0600 (kassick)"
 * $Id$
 * Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
 */
/*
 * ===========================================================================
 *
 *       Filename:  uniq_ids_q.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  13-11-2011 20:56:10 BRST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include "src/libRastro/list.h"
#include "src/libRastro/uniq_ids_q.h"

// get config struct
#include "server-config.h"
#include "src/server/pvfs2-server.h"

#include <rastro.h>

unique_ids_queue_t unique_io_ids;
int rst_server_id;

void init_unique_op_id(unique_ids_queue_t * q)
{
  INIT_LIST_HEAD(& (q->free_list) );
  pthread_mutex_init(& (q->mutex), NULL );
  q->max = 0;
}

int get_unique_op_id(unique_ids_queue_t * q, rst_buffer_t ** buf)
{
  int ret;

  pthread_mutex_lock(& (q->mutex)  );

  if (list_empty( &( q->free_list) ))
  {
    // new op, get a new rastro buffer and initialize it, but DO NOT SET
    // THE DESTROY CALLBACK

    ret = q->max ++;
    *buf = (rst_buffer_t *) malloc(sizeof(rst_buffer_t));
    // create file and prepare everything, but do not run set_specific --
    // and just in case, do not allow any stray destructor do release the
    // buffer
    (*buf)->do_destroy = 0;
    rst_init_ptr_set_key(*buf, rst_server_id, 100 + ret , 0);

  } else {
    // use the first op available on the list
    unique_id_entry_t *tmp = list_entry(q->free_list.next,
                                      unique_id_entry_t, FromAndToList );
    *buf = tmp->rst_ptr;
    ret = tmp->id;
    list_del( &(tmp->FromAndToList) );
    free(tmp);
  }

  pthread_mutex_unlock( &(q->mutex) );
  return ret;
}


void release_unique_op_id(unique_ids_queue_t * q, int id, rst_buffer_t * buf)
{

  // create a list entry to store this id
  unique_id_entry_t *tmp = malloc(sizeof(unique_id_entry_t));
  tmp->id = id;
  tmp->rst_ptr = buf;

  pthread_mutex_lock(& (q->mutex)  );

  list_add_tail( &(tmp->FromAndToList), &(q->free_list) );

  pthread_mutex_unlock( &(q->mutex) );
}

