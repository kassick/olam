/* C source code
 * File: "/home/kassick/Work/olam/olamfs-trace-converter/src/olam_trace2paje.c"
 * Created: "Ter, 31 Mai 2011 11:11:38 -0300 (kassick)"
 * Updated: "Ter, 07 Jun 2011 18:41:26 -0300 (kassick)"
 * $Id$
 * Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
 */
/*
 * ===========================================================================
 *
 *       Filename:  olam_trace2paje.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  31-05-2011 11:11:38 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */

#include <stdlib.h>
#include <getopt.h>
#include <stdio.h>
#include <search.h>
#include <rastro.h>
#include <assert.h>

#include <aky.h>
#include <pvfs_events.h>
#include <olam_events.h>

#define FALSE (0)
#define TRUE (!FALSE)

#define OPTSTR "mho"
#define PROGNAME "olam_trace2paje"


void usage() {
  printf("Usage: %s <options> <file1.rst> <file2.rst> ...\n",PROGNAME);
  printf(" -h                 displays this message\n");
  printf(" -m                 adds MPI messages to the paje file\n");
  printf(" -o                 disables offset in constants (for olam rastro v>=11)");
  printf(" -somethingels      does some thing else\n");
}

int main(int argc, char** argv)
{
  int opt;
  char mpi_messages = FALSE;
  int const_offset = OLAM_EVT_BASE;
  
  rst_file_t data;
  rst_event_t event;
  int i;



  while ((opt = getopt(argc,argv, OPTSTR) != -1 ))
  {
    switch(opt) {
      case 'm':
        mpi_messages = TRUE;
        break;

      case 'h':
        usage();
        exit(0);

      case 'o':
        const_offset = 0;
        break;

      default:
        printf("Unknown option -%c\n",(char)opt);
        usage();
        exit(1);
    }
  }

  if ((argc - optind) == 0)
  {
    usage();
    exit(1);
  }


  // ---
  // Done with all the bureocracies


  hcreate(1000000);

  for (i = optind; i < argc; i++) {
    int ret = rst_open_file(argv[i], &data, NULL, 100000);
    if (ret == -1) {
      printf("%s: trace %s could not be opened\n", argv[0], argv[i]);
      return 1;
    }
  }

  name_init();
  paje_header();
  paje_hierarchy();
  pajeCreateContainer(0,"OLAM","OLAM","0","OLAM");
  pajeCreateContainer(0,"FS","FS","0","FS");

  //printf("Hello!\n");

  while (rst_decode_event(&data, &event)) {
    char mpi_process[100];
    char mpi_name[100];
    char value[100];
    char state[100];
    double timestamp;
    int type;
    char key[AKY_DEFAULT_STR_SIZE];
    char *locfn,*locfn_access;

    if (event.id2 >= PVFS_VERSION_BASE)
    {
      // this is a PVFS event
      // Identified by server id and by file name
      snprintf(mpi_process, 100, "pvfs_%ld",event.id1);
    } else { 
      // this is olam event -- identified by the process number
      snprintf(mpi_process, 100, "rank%ld", event.id1);
    }
    
    snprintf(value, 100, "%s", name_get(event.type));
    timestamp = (double) event.timestamp / 1000000;


    switch (event.type) {
      // Use Aky events as well
      case MPI_INIT:
      case PVFS_EVT_INIT:
        pajeCreateContainer(timestamp, mpi_process,
                            "PROCESS", "0", mpi_process);
        break;
      case OLAM_EVT_HDF5_OPEN_IN:
        // Creates the container /filename/ 
        // Creates arrow from rank%d to filename
        // Creates an state "opened"
        // open is iiss, myrank, thread_id, locfn,access
        assert(event.ct.n_string == 2);
        locfn = event.v_string[0];
        locfn_access = event.v_string[1];
      
        // Type here is FILE -- Hope paje does not go crazy... ;)
        pajeCreateContainer(timestamp, locfn,
                            "FILE", "0", locfn);
        
        snprintf(state,100, "STATE_%s",value);
        pajePushState(timestamp, locfn, state, value);
        snprintf(state,100, "STATE_%s","opened");
        pajePushState(timestamp, locfn, state, "opened");

        // id is "open W ## myrank ## /file/name"
        snprintf(key,AKY_DEFAULT_STR_SIZE,"open %s#%d#%s",
            locfn_access, event.id1, locfn);
        pajeStartLink(timestamp, "0", "LINK", mpi_process, "PTP", key);
        pajeEndLink  (timestamp, "0", "LINK", locfn      , "PTP", key);
        
        break;


      case OLAM_EVT_HDF5_CLOSE_IN:
      case OLAM_EVT_HDF5_CREATE_IN:
        // Creates arrow from rank%d to filename
        // Creates an state with event

        assert(event.ct.n_string == 2);
        locfn = event.v_string[0];
        //locfn_access = event.v_string[1];
       
        snprintf(state,100, "STATE_%s",value);
        pajePushState(timestamp, locfn, state, value);

        // id is "open W ## myrank ## /file/name"
        snprintf(key,AKY_DEFAULT_STR_SIZE,"%s#%d#%s",
            mpi_name, event.id1, locfn);
        pajeStartLink(timestamp, "0", "LINK", mpi_process, "PTP", key);
        pajeEndLink  (timestamp, "0", "LINK", locfn      , "PTP", key);
        
        break;

      case OLAM_EVT_HDF5_OPEN_OUT:
      case OLAM_EVT_HDF5_CREATE_OUT:
        // Closes an state "open"
        // open is iiss, myrank, thread_id, locfn,access
        assert(event.ct.n_string == 2);
        locfn = event.v_string[0];
        locfn_access = event.v_string[1];
      
        snprintf(state,100, "STATE_%s",value);
        pajePopState(timestamp, locfn, state);

        break;


      case OLAM_EVT_HDF5_CLOSE_OUT:
        // Destroys the container /filename/ 
        // Creates arrow from rank%d to filename
        // closes an state "opened"
        // open is iiss, myrank, thread_id, locfn,access
        assert(event.ct.n_string == 2);
        locfn = event.v_string[0];
      
        // Type here is FILE -- Hope paje does not go crazy... ;)
        pajeCreateContainer(timestamp, locfn,
                            "FILE", "0", locfn);
        
        snprintf(state,100, "STATE_%s",value);
        pajePushState(timestamp, locfn, state, value);
        snprintf(state,100, "STATE_%s","opened");
        pajePushState(timestamp, locfn, state, "opened");

        // id is "open W ## myrank ## /file/name"
        snprintf(key,AKY_DEFAULT_STR_SIZE,"open %s#%d#%s",
            locfn_access, event.id1, locfn);
        pajeStartLink(timestamp, "0", "LINK", mpi_process, "PTP", key);
        pajeEndLink  (timestamp, "0", "LINK", locfn      , "PTP", key);
        
        break;

      case MPI_COMM_SPAWN_IN:
      case MPI_COMM_GET_NAME_IN:
      case MPI_COMM_SET_NAME_IN:
      case MPI_REDUCE_IN:
      case MPI_ALLREDUCE_IN:
      case MPI_REDUCE_SCATTER_IN:
      case MPI_ALLGATHER_IN:
      case MPI_ALLGATHERV_IN:
      case MPI_SCATTER_IN:
      case MPI_SCATTERV_IN:
      case MPI_WAIT_IN:
      case MPI_IRECV_IN:
      case MPI_ISEND_IN:
      case MPI_RECV_IN:
      case MPI_SEND_IN:
      case MPI_BCAST_IN:
      case MPI_BARRIER_IN:
      case MPI_GATHER_IN:
      case MPI_GATHERV_IN:
      case MPI_ALLTOALL_IN:
      case MPI_ALLTOALLV_IN:
      case MPI_OP_CREATE_IN:
      case MPI_OP_FREE_IN:
      case MPI_SCAN_IN:
      case MPI_ATTR_DELETE_IN:
      case MPI_ATTR_GET_IN:
      case MPI_ATTR_PUT_IN:
      case MPI_COMM_COMPARE_IN:
      case MPI_COMM_CREATE_IN:
      case MPI_COMM_DUP_IN:
      case MPI_COMM_FREE_IN:
      case MPI_COMM_GROUP_IN:
      case MPI_COMM_RANK_IN:
      case MPI_COMM_REMOTE_GROUP_IN:
      case MPI_COMM_REMOTE_SIZE_IN:
      case MPI_COMM_SIZE_IN:
      case MPI_COMM_SPLIT_IN:
      case MPI_COMM_TEST_INTER_IN:
      case MPI_GROUP_COMPARE_IN:
      case MPI_GROUP_DIFFERENCE_IN:
      case MPI_GROUP_EXCL_IN:
      case MPI_GROUP_FREE_IN:
      case MPI_GROUP_INCL_IN:
      case MPI_GROUP_INTERSECTION_IN:
      case MPI_GROUP_RANK_IN:
      case MPI_GROUP_RANGE_EXCL_IN:
      case MPI_GROUP_RANGE_INCL_IN:
      case MPI_GROUP_SIZE_IN:
      case MPI_GROUP_TRANSLATE_RANKS_IN:
      case MPI_GROUP_UNION_IN:
      case MPI_INTERCOMM_CREATE_IN:
      case MPI_INTERCOMM_MERGE_IN:
      case MPI_KEYVAL_CREATE_IN:
      case MPI_KEYVAL_FREE_IN:
      case MPI_ABORT_IN:
      case MPI_ERROR_CLASS_IN:
      case MPI_ERRHANDLER_CREATE_IN:
      case MPI_ERRHANDLER_FREE_IN:
      case MPI_ERRHANDLER_GET_IN:
      case MPI_ERROR_STRING_IN:
      case MPI_ERRHANDLER_SET_IN:
      case MPI_GET_PROCESSOR_NAME_IN:
      case MPI_INITIALIZED_IN:
      case MPI_WTICK_IN:
      case MPI_WTIME_IN:
      case MPI_ADDRESS_IN:
      case MPI_BSEND_IN:
      case MPI_BSEND_INIT_IN:
      case MPI_BUFFER_ATTACH_IN:
      case MPI_BUFFER_DETACH_IN:
      case MPI_CANCEL_IN:
      case MPI_REQUEST_FREE_IN:
      case MPI_RECV_INIT_IN:
      case MPI_SEND_INIT_IN:
      case MPI_GET_ELEMENTS_IN:
      case MPI_GET_COUNT_IN:
      case MPI_IBSEND_IN:
      case MPI_IPROBE_IN:
      case MPI_IRSEND_IN:
      case MPI_ISSEND_IN:
      case MPI_PACK_IN:
      case MPI_PACK_SIZE_IN:
      case MPI_PROBE_IN:
      case MPI_RSEND_IN:
      case MPI_RSEND_INIT_IN:
      case MPI_SENDRECV_IN:
      case MPI_SENDRECV_REPLACE_IN:
      case MPI_SSEND_IN:
      case MPI_SSEND_INIT_IN:
      case MPI_START_IN:
      case MPI_STARTALL_IN:
      case MPI_TEST_IN:
      case MPI_TESTALL_IN:
      case MPI_TESTANY_IN:
      case MPI_TEST_CANCELLED_IN:
      case MPI_TESTSOME_IN:
      case MPI_TYPE_COMMIT_IN:
      case MPI_TYPE_CONTIGUOUS_IN:
      case MPI_TYPE_EXTENT_IN:
      case MPI_TYPE_FREE_IN:
      case MPI_TYPE_HINDEXED_IN:
      case MPI_TYPE_HVECTOR_IN:
      case MPI_TYPE_INDEXED_IN:
      case MPI_TYPE_LB_IN:
      case MPI_TYPE_SIZE_IN:
      case MPI_TYPE_STRUCT_IN:
      case MPI_TYPE_UB_IN:
      case MPI_TYPE_VECTOR_IN:
      case MPI_UNPACK_IN:
      case MPI_WAITALL_IN:
      case MPI_WAITANY_IN:
      case MPI_WAITSOME_IN:
      case MPI_CART_COORDS_IN:
      case MPI_CART_CREATE_IN:
      case MPI_CART_GET_IN:
      case MPI_CART_MAP_IN:
      case MPI_CART_SHIFT_IN:
      case MPI_CARTDIM_GET_IN:
      case MPI_DIMS_CREATE_IN:
      case MPI_GRAPH_CREATE_IN:
      case MPI_GRAPH_GET_IN:
      case MPI_GRAPH_MAP_IN:
      case MPI_GRAPH_NEIGHBORS_IN:
      case MPI_GRAPH_NEIGHBORS_COUNT_IN:
      case MPI_GRAPHDIMS_GET_IN:
      case MPI_TOPO_TEST_IN:
      case MPI_RECV_IDLE_IN:
      case MPI_CART_RANK_IN:
      case MPI_CART_SUB_IN:
      case MPI_FINALIZE_IN:
        pajePushState(timestamp, mpi_process, "STATE", value);
        break;
      case MPI_COMM_SPAWN_OUT:
      case MPI_COMM_GET_NAME_OUT:
      case MPI_COMM_SET_NAME_OUT:
      case MPI_REDUCE_OUT:
      case MPI_ALLREDUCE_OUT:
      case MPI_REDUCE_SCATTER_OUT:
      case MPI_ALLGATHER_OUT:
      case MPI_ALLGATHERV_OUT:
      case MPI_SCATTER_OUT:
      case MPI_SCATTERV_OUT:
      case MPI_WAIT_OUT:
      case MPI_IRECV_OUT:
      case MPI_ISEND_OUT:
      case MPI_RECV_OUT:
      case MPI_SEND_OUT:
      case MPI_BCAST_OUT:
      case MPI_BARRIER_OUT:
      case MPI_GATHER_OUT:
      case MPI_GATHERV_OUT:
      case MPI_ALLTOALL_OUT:
      case MPI_ALLTOALLV_OUT:
      case MPI_OP_CREATE_OUT:
      case MPI_OP_FREE_OUT:
      case MPI_SCAN_OUT:
      case MPI_ATTR_DELETE_OUT:
      case MPI_ATTR_GET_OUT:
      case MPI_ATTR_PUT_OUT:
      case MPI_COMM_COMPARE_OUT:
      case MPI_COMM_CREATE_OUT:
      case MPI_COMM_DUP_OUT:
      case MPI_COMM_FREE_OUT:
      case MPI_COMM_GROUP_OUT:
      case MPI_COMM_RANK_OUT:
      case MPI_COMM_REMOTE_GROUP_OUT:
      case MPI_COMM_REMOTE_SIZE_OUT:
      case MPI_COMM_SIZE_OUT:
      case MPI_COMM_SPLIT_OUT:
      case MPI_COMM_TEST_INTER_OUT:
      case MPI_GROUP_COMPARE_OUT:
      case MPI_GROUP_DIFFERENCE_OUT:
      case MPI_GROUP_EXCL_OUT:
      case MPI_GROUP_FREE_OUT:
      case MPI_GROUP_INCL_OUT:
      case MPI_GROUP_INTERSECTION_OUT:
      case MPI_GROUP_RANK_OUT:
      case MPI_GROUP_RANGE_EXCL_OUT:
      case MPI_GROUP_RANGE_INCL_OUT:
      case MPI_GROUP_SIZE_OUT:
      case MPI_GROUP_TRANSLATE_RANKS_OUT:
      case MPI_GROUP_UNION_OUT:
      case MPI_INTERCOMM_CREATE_OUT:
      case MPI_INTERCOMM_MERGE_OUT:
      case MPI_KEYVAL_CREATE_OUT:
      case MPI_KEYVAL_FREE_OUT:
      case MPI_ABORT_OUT:
      case MPI_ERROR_CLASS_OUT:
      case MPI_ERRHANDLER_CREATE_OUT:
      case MPI_ERRHANDLER_FREE_OUT:
      case MPI_ERRHANDLER_GET_OUT:
      case MPI_ERROR_STRING_OUT:
      case MPI_ERRHANDLER_SET_OUT:
      case MPI_GET_PROCESSOR_NAME_OUT:
      case MPI_INITIALIZED_OUT:
      case MPI_WTICK_OUT:
      case MPI_WTIME_OUT:
      case MPI_ADDRESS_OUT:
      case MPI_BSEND_OUT:
      case MPI_BSEND_INIT_OUT:
      case MPI_BUFFER_ATTACH_OUT:
      case MPI_BUFFER_DETACH_OUT:
      case MPI_CANCEL_OUT:
      case MPI_REQUEST_FREE_OUT:
      case MPI_RECV_INIT_OUT:
      case MPI_SEND_INIT_OUT:
      case MPI_GET_ELEMENTS_OUT:
      case MPI_GET_COUNT_OUT:
      case MPI_IBSEND_OUT:
      case MPI_IPROBE_OUT:
      case MPI_IRSEND_OUT:
      case MPI_ISSEND_OUT:
      case MPI_PACK_OUT:
      case MPI_PACK_SIZE_OUT:
      case MPI_PROBE_OUT:
      case MPI_RSEND_OUT:
      case MPI_RSEND_INIT_OUT:
      case MPI_SENDRECV_OUT:
      case MPI_SENDRECV_REPLACE_OUT:
      case MPI_SSEND_OUT:
      case MPI_SSEND_INIT_OUT:
      case MPI_START_OUT:
      case MPI_STARTALL_OUT:
      case MPI_TEST_OUT:
      case MPI_TESTALL_OUT:
      case MPI_TESTANY_OUT:
      case MPI_TEST_CANCELLED_OUT:
      case MPI_TESTSOME_OUT:
      case MPI_TYPE_COMMIT_OUT:
      case MPI_TYPE_CONTIGUOUS_OUT:
      case MPI_TYPE_EXTENT_OUT:
      case MPI_TYPE_FREE_OUT:
      case MPI_TYPE_HINDEXED_OUT:
      case MPI_TYPE_HVECTOR_OUT:
      case MPI_TYPE_INDEXED_OUT:
      case MPI_TYPE_LB_OUT:
      case MPI_TYPE_SIZE_OUT:
      case MPI_TYPE_STRUCT_OUT:
      case MPI_TYPE_UB_OUT:
      case MPI_TYPE_VECTOR_OUT:
      case MPI_UNPACK_OUT:
      case MPI_WAITALL_OUT:
      case MPI_WAITANY_OUT:
      case MPI_WAITSOME_OUT:
      case MPI_CART_COORDS_OUT:
      case MPI_CART_CREATE_OUT:
      case MPI_CART_GET_OUT:
      case MPI_CART_MAP_OUT:
      case MPI_CART_SHIFT_OUT:
      case MPI_CARTDIM_GET_OUT:
      case MPI_DIMS_CREATE_OUT:
      case MPI_GRAPH_CREATE_OUT:
      case MPI_GRAPH_GET_OUT:
      case MPI_GRAPH_MAP_OUT:
      case MPI_GRAPH_NEIGHBORS_OUT:
      case MPI_GRAPH_NEIGHBORS_COUNT_OUT:
      case MPI_GRAPHDIMS_GET_OUT:
      case MPI_TOPO_TEST_OUT:
      case MPI_RECV_IDLE_OUT:
      case MPI_CART_RANK_OUT:
      case MPI_CART_SUB_OUT:
        pajePopState(timestamp, mpi_process, "STATE");
        break;
        locfn = event.v_string[0];
        pajePopState(timestamp, locfn , "STATE");
      case MPI_FINALIZE_OUT:
        pajePopState(timestamp, mpi_process, "STATE");
        pajeDestroyContainer(timestamp, "PROCESS", mpi_process);
        break;
      

#if 0
    case AKY_PTP_SEND:
      {
        char key[AKY_DEFAULT_STR_SIZE];
        aky_put_key("n", event.id1, event.v_uint32[0], key,
                    AKY_DEFAULT_STR_SIZE);
        pajeStartLink(timestamp, "0", "LINK", mpi_process, "PTP", key);
      }
      break;
    case AKY_PTP_RECV:
      {
        char key[AKY_DEFAULT_STR_SIZE];
        aky_get_key("n", event.v_uint32[0], event.id1, key,
                    AKY_DEFAULT_STR_SIZE);
        pajeEndLink(timestamp, "0", "LINK", mpi_process, "PTP", key);
      }
      break;
#endif
    }
  }

  rst_close_file(&data);
  hdestroy();
  return 0;

}
