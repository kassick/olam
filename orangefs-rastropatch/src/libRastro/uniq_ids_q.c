/* C source code
 * File: "/home/kassick/Work/orangefs-rastropatch/src/libRastro/uniq_ids_q.c"
 * Created: "Dom, 13 Nov 2011 20:56:10 -0200 (kassick)"
 * Updated: "Dom, 13 Nov 2011 23:03:21 -0200 (kassick)"
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

unique_ids_queue_t unique_io_ids;

void init_unique_op_id(unique_ids_queue_t * q)
{
  INIT_LIST_HEAD(& (q->free_list) );
  pthread_mutex_init(& (q->mutex), NULL );
  q->max = 0;
}

int get_unique_op_id(unique_ids_queue_t * q)
{
  int ret;
  pthread_mutex_lock(& (q->mutex)  );

  if (list_empty( &( q->free_list) ))
  {
    ret = q->max++;
  } else {
    unique_id_entry_t *tmp = list_entry(q->free_list.next,
                                      unique_id_entry_t, FromAndToList );
    ret = tmp->id;
    list_del( &(tmp->FromAndToList) );
    free(tmp);
  }

  pthread_mutex_unlock( &(q->mutex) );
  return ret;
}
void release_unique_op_id(unique_ids_queue_t * q, int id)
{
  pthread_mutex_lock(& (q->mutex)  );

  if (id == q->max-1) {
    // if possible, just decrease the max
    q->max--;
  } else {
    // create a list entry to store this id
    unique_id_entry_t *tmp = malloc(sizeof(unique_id_entry_t));
    tmp->id = id;

    list_add_tail( &(tmp->FromAndToList), &(q->free_list) );
  }
  pthread_mutex_unlock( &(q->mutex) );
}

