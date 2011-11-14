// C++ source code
// File: "/home/kassick/Work/orangefs-rastropatch/src/libRastro/uniq_ids_q.h"
// Created: "Dom, 13 Nov 2011 20:56:48 -0200 (kassick)"
// Updated: "Dom, 13 Nov 2011 23:00:44 -0200 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  uniq_ids_q.h
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  13-11-2011 20:56:48 BRST
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */


#include "src/libRastro/list.h"
#include <pthread.h>

typedef struct { 
  struct list_head free_list;
  int max;
  pthread_mutex_t mutex;

} unique_ids_queue_t;

typedef struct _unique_id_entry_t {
    struct list_head FromAndToList;
    int id;
} unique_id_entry_t;



void init_unique_op_id(unique_ids_queue_t * q);
int get_unique_op_id(unique_ids_queue_t * q);
void release_unique_op_id(unique_ids_queue_t * q, int id);


extern unique_ids_queue_t unique_io_ids;

