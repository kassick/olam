// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/logging.cc"
// Created: "Qui, 24 Nov 2011 21:24:24 -0600 (kassick)"
// Updated: "Qui, 24 Nov 2011 21:40:32 -0600 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  logging.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  24-11-2011 21:24:24 CST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include "logging.hh"


namespace logging {
const char *  log_preceed(unsigned int level) {
        switch (level) {
          case FLAG_INFO:
            return "Information";
          case FLAG_WARN:
            return "Warning";
          case FLAG_ERROR:
            return "Error";
          case FLAG_FATAL:
            return "Fatal";
        }

        return "unknown";
  }

}
