// C++ source code
// File: "/home/kassick/Work/olam/olamfs-trace-converter/include/pvfs_events.h"
// Created: "Sex, 03 Jun 2011 14:42:25 -0300 (kassick)"
// Updated: "Ter, 07 Jun 2011 18:28:17 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  pvfs_events.h
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  03-06-2011 14:42:25 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */

#ifndef _PVFS_EVENTS_H_
#define _PVFS_EVENTS_H_

#define PVFS_VERSION_BASE (20)
#define PVFS_EVT_BASE (4500)

// IN / OUT Declaration
#define PVFS_EVT_IO(name) \
 PVFS_EVT_ ## name ## _IN,\
 PVFS_EVT_ ## name ## _OUT

#define PVFS_EVT_N(name) \
  PVFS_EVT_ ## name

#define PVFS_EVT_STR(id) \
  {PVFS_EVT_ ## id ## _IN, "PVFS_" #id}, \
  {PVFS_EVT_ ## id ## _OUT, "PVFS_" #id "V"}

#define PVFS_EVT_STR_N(id) \
  {PVFS_EVT_ ## id  , "PVFS_" #id}



// Add IDs here
typedef enum _pvfs_evts_t {
  PVFS_EVT_INIT = PVFS_EVT_BASE,
  PVFS_EVT_IO(CREATE),
  PVFS_EVT_IO(OPEN),
  PVFS_EVT_IO(WRITE),
  PVFS_EVT_IO(READ),
} pvfs_evts_t;




typedef struct _pvfs_evt_names_t {
  pvfs_evts_t id;
  char * name;
} pvfs_evt_names_t;

// Add names to list here
static pvfs_evt_names_t  pvfs_evt_names[] = { 
  PVFS_EVT_STR_N(INIT),
  PVFS_EVT_STR(OPEN),
  PVFS_EVT_STR(CREATE),
  PVFS_EVT_STR(WRITE),
  PVFS_EVT_STR(READ),
  {0, NULL},
};


#endif
