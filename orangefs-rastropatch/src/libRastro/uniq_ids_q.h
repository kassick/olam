// C++ source code
// File: "/home/kassick/Work/olam/orangefs-rastropatch/src/libRastro/uniq_ids_q.h"
// Created: "Dom, 13 Nov 2011 20:56:48 -0200 (kassick)"
// Updated: "Sáb, 26 Nov 2011 00:34:16 -0600 (kassick)"
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


#ifndef __UNIQ_TYPES_H__
#define __UNIQ_TYPES_H__

#include "src/libRastro/list.h"
#include <pthread.h>
#include <rastro.h>

typedef struct { 
  struct list_head free_list;
  int max;
  pthread_mutex_t mutex;

} unique_ids_queue_t;

typedef struct _unique_id_entry_t {
    struct list_head FromAndToList;
    int id;

    // a per-op rastro buffer because ot mtfkng threads
    rst_buffer_t *rst_ptr;

} unique_id_entry_t;



void init_unique_op_id(unique_ids_queue_t * q);
int get_unique_op_id(unique_ids_queue_t * q, rst_buffer_t ** buf);
void release_unique_op_id(unique_ids_queue_t * q, int id, rst_buffer_t * buf);


// this should be elseqwhere to make this more generic...
extern unique_ids_queue_t unique_io_ids;
extern int rst_server_id;


#endif
