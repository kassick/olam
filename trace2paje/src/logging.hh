// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/logging.hh"
// Created: "Ter, 22 Nov 2011 16:03:50 -0600 (kassick)"
// Updated: "Qui, 24 Nov 2011 21:28:45 -0600 (kassick)"
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
 *        Created:  22-11-2011 16:03:50 CST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#ifndef __LOGGING_HH_INCLUDED__
#define __LOGGING_HH_INCLUDED__

#include <fstream>
#include <iostream>
      
#define  FLAG_NONE  ( 2 << 0)
#define  FLAG_INFO  ( 2 << 1)
#define  FLAG_WARN  ( 2 << 2)
#define  FLAG_ERROR ( 2 << 3)
#define  FLAG_FATAL ( 2 << 4)


#define LOGIFLEVEL(logger,alevel)  \
  if (!(logger.level & alevel)) ;\
  else logger.fout  << "[" << logging::log_preceed(alevel) << "] "

#define WARN(logger)  LOGIFLEVEL(logger,FLAG_WARN)
#define INFO(logger)  LOGIFLEVEL(logger,FLAG_INFO)
#define ERROR(logger)  LOGIFLEVEL(logger,FLAG_ERROR)
#define FATAL(logger)  LOGIFLEVEL(logger,FLAG_FATAL)

using namespace std;
namespace logging {


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Return the string to be shown for a debug level
 *  Description:  
 * =====================================================================================
 */

const char *  log_preceed(unsigned int level);


  class Logger {
    public:

      ofstream fout;

      unsigned int level;


      Logger() {
        fout.open("/dev/stderr");
      }

      Logger (unsigned int _level)
      {
        level = _level;
      }

      Logger (unsigned int _level, const string & logfile)
      {
        level = _level;
        fout.open(logfile);

        if (!fout.is_open())
        {
          cerr << "Could not open specified logfile, using stderr" << endl;
        }
      }

  };
}



extern logging::Logger logger;

#endif
