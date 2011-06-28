// C++ source code
// File: "/home/kassick/Work/olam/olamfs-trace-converter/include/pvfs_events.h"
// Created: "Sex, 03 Jun 2011 14:42:25 -0300 (kassick)"
// Updated: "Sex, 17 Jun 2011 18:11:22 -0300 (kassick)"
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

#include <stdlib.h>
#include <name_types.h>
#include <macros_f.h>

#define PVFS_VERSION_BASE (20)
#define PVFS_EVT_BASE (4500)
#define PVFS_EVT_MAX  (4999)
#define EVT_BASE_NAME PVFS
#include <evt_template.h>

// Add IDs here
typedef enum _pvfs_evts_t {
  PVFS_INIT = PVFS_EVT_BASE,
  EVT_IO(CREATE),
  EVT_IO(OPEN),
  EVT_IO(WRITE),
  EVT_IO(READ),
} pvfs_evts_t;



// Add names to list here
static evt_name_t  pvfs_evt_names[] = { 
  EVT_NAME_ENTRY_N(INIT),
  EVT_NAME_ENTRY(OPEN),
  EVT_NAME_ENTRY(CREATE),
  EVT_NAME_ENTRY(WRITE),
  EVT_NAME_ENTRY(READ),
  NULL_EVT,
};


#undef EVT_BASE_NAME
#endif
