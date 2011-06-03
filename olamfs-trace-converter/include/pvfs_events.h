// C++ source code
// File: "/home/kassick/Work/olam/olamfs-trace-converter/include/pvfs_events.h"
// Created: "Sex, 03 Jun 2011 14:42:25 -0300 (kassick)"
// Updated: "Sex, 03 Jun 2011 15:31:42 -0300 (kassick)"
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

#define _PVFS_CUR_EVT (PVFS_EVT_BASE)

#define PVFS_EVT_IO(name) \\
#define PVFS_EVT_ ## name ## _IN (_PVFS_CUR_EVT) \\
#define PVFS_EVT_ ## name ## _OUT (_PVFS_CUR_EVT + 1) \\
#define _PVFS_CUR_EVT (_PVFS_CUR_EVT + 2 )

#define PVFS_EVT_N(name) \\
#define PVFS_EVT_ ## name ## _IN (_PVFS_CUR_EVT) \\
#define _PVFS_CUR_EVT (_PVFS_CUR_EVT + 1 )

PVFS_EVT_N(INIT)
PVFS_EVT_IO(WRITE)
PVFS_EVT_IO(READ)


#endif
