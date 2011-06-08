// C++ source code
// File: "/home/kassick/Work/olam/olamfs-trace-converter/include/macros_f.h"
// Created: "Qua, 08 Jun 2011 18:58:37 -0300 (kassick)"
// Updated: "Qua, 08 Jun 2011 19:07:46 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  macros_f.h
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  08-06-2011 18:58:37 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */

#ifndef MACROS_F_H
#define MACROS_F_H

#define xstr(s) str(s)
#define str(s) #s
#define nop(s) s
#define concat(s1,s2) s1 ## s2
#define CAT1(X,Y) X##_##Y
#define CAT3(X,Y,Z) X##_##Y##_##Z
#define CAT4(X,Y,Z,W) X##_##Y##_##Z##_##W
#define TEMPLATE(X,Y) CAT1(X,Y)
#define TEMPLATE3(X,Y,Z) CAT3(X,Y,Z)
#define TEMPLATE4(X,Y,Z,W) CAT4(X,Y,Z,W)


#endif
