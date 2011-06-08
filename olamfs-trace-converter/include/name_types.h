// C++ source code
// File: "/home/kassick/Work/olam/olamfs-trace-converter/include/name_types.h"
// Created: "Qua, 08 Jun 2011 18:03:19 -0300 (kassick)"
// Updated: "Qua, 08 Jun 2011 19:03:23 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  name_types.h
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  08-06-2011 18:03:19 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */


#ifndef _NAME_TYPES_H_
#define _NAME_TYPES_H_

#define NULL_EVT {0, 0, NULL, NULL, NULL, EVT}
// Types can be defined only once, but the macros need to be ?
typedef enum _evt_type_t {
  IN = 0,
  OUT,
  EVT
} evt_type_t;

typedef struct _evt_names_t {
  int id,start_id;
  char * name, *short_name,*start_name;
  int type;
} evt_name_t;

#endif
