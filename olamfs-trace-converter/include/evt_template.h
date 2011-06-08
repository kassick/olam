// C++ source code
// File: "/home/kassick/Work/olam/olamfs-trace-converter/include/evt_template.h"
// Created: "Qua, 08 Jun 2011 19:00:14 -0300 (kassick)"
// Updated: "Qua, 08 Jun 2011 19:10:46 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  evt_template.h
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  08-06-2011 19:00:14 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */

#define EVT_IO(name) \
 TEMPLATE4(nop(EVT_BASE_NAME),EVT,name,IN), \
 TEMPLATE4(nop(EVT_BASE_NAME),EVT,name,OUT)

#define EVT_N(name) \
  TEMPLATE3(nop(EVT_BASE_NAME),EVT,name)

#define EVT_NAME_ENTRY(id) \
  {TEMPLATE4(nop(EVT_BASE_NAME),EVT,id,IN) , TEMPLATE4(nop(EVT_BASE_NAME),EVT,id,IN), xstr(EVT_BASE_NAME) "_" #id     , #id, xstr(EVT_BASE_NAME) "_" #id, IN}, \
  {TEMPLATE4(nop(EVT_BASE_NAME),EVT,id,OUT), TEMPLATE4(nop(EVT_BASE_NAME),EVT,id,IN), xstr(EVT_BASE_NAME) "_" #id "V" , #id, xstr(EVT_BASE_NAME) "_" #id, OUT}

#define EVT_NAME_ENTRY_N(id) \
  {TEMPLATE3(nop(EVT_BASE_NAME),EVT,id)    , TEMPLATE3(nop(EVT_BASE_NAME),EVT,id)   , xstr(EVT_BASE_NAME) "_" #id     , #id, xstr(EVT_BASE_NAME) "_" #id, EVT} 
  //{OLAM_EVT_ ## id        , OLAM_EVT_ ## id        ,"OLAM_" #id     , #id, "OLAM_" #id, EVT}


