// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/rastro_helper.hh"
// Created: "Sex, 16 Set 2011 18:23:38 -0300 (kassick)"
// Updated: "Sex, 16 Set 2011 18:25:12 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  rastro_helper.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  16-09-2011 18:23:38 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */


#ifndef __RASTRO_HELPER_H__
#define __RASTRO_HELPER_H__

typedef enum {
    c, w, i, l, f, d, s, invalid,
} rastro_basic_types_t;

typedef union {
        char c;
        uint16_t w;
        uint32_t i;
        uint64_t l;
        float f;
        double d;
        char *s;
} rastro_basic_val_t;




#endif
