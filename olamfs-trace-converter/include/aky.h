/*
    This file is part of Akypuera

    Akypuera is free software: you can redistribute it and/or modify
    it under the terms of the GNU Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Akypuera is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Public License for more details.

    You should have received a copy of the GNU Public License
    along with Akypuera. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef __AKY_H
#define __AKY_H

#define MPI_INIT 1000
#define MPI_FINALIZE_IN 1002
#define MPI_FINALIZE_OUT 1003

/* MPI-2 */
#define MPI_COMM_SPAWN_IN 2000
#define MPI_COMM_SPAWN_OUT 2001

/* Additional Communicators */
#define MPI_COMM_GET_NAME_IN 3000
#define MPI_COMM_GET_NAME_OUT 3001
#define MPI_COMM_SET_NAME_IN 3002
#define MPI_COMM_SET_NAME_OUT 3003

#define MAX_AKY_STATE_NAMES 5002

/* used */
#define MPI_REDUCE_IN 25
#define MPI_REDUCE_OUT 26
#define MPI_ALLREDUCE_IN 5
#define MPI_ALLREDUCE_OUT 6
#define MPI_REDUCE_SCATTER_IN 23
#define MPI_REDUCE_SCATTER_OUT 24
#define MPI_ALLGATHER_IN 1
#define MPI_ALLGATHER_OUT 2
#define MPI_ALLGATHERV_IN 3
#define MPI_ALLGATHERV_OUT 4
#define MPI_SCATTER_IN 29
#define MPI_SCATTER_OUT 30
#define MPI_SCATTERV_IN 31
#define MPI_SCATTERV_OUT 32
#define MPI_WAIT_IN 211
#define MPI_WAIT_OUT 212
#define MPI_IRECV_IN 141
#define MPI_IRECV_OUT 142
#define MPI_ISEND_IN 145
#define MPI_ISEND_OUT 146
#define MPI_RECV_IN 155
#define MPI_RECV_OUT 156
#define MPI_SEND_IN 161
#define MPI_SEND_OUT 162
#define MPI_BCAST_IN 13
#define MPI_BCAST_OUT 14
#define MPI_BARRIER_IN 11
#define MPI_BARRIER_OUT 12
#define MPI_GATHER_IN 15
#define MPI_GATHER_OUT 16
#define MPI_GATHERV_IN 17
#define MPI_GATHERV_OUT 18
#define MPI_ALLTOALL_IN 7
#define MPI_ALLTOALL_OUT 8
#define MPI_ALLTOALLV_IN 9
#define MPI_ALLTOALLV_OUT 10

/* not used */
#define MPI_OP_CREATE_IN 19
#define MPI_OP_CREATE_OUT 20
#define MPI_OP_FREE_IN 21
#define MPI_OP_FREE_OUT 22
#define MPI_SCAN_IN 27
#define MPI_SCAN_OUT 28
#define MPI_ATTR_DELETE_IN 33
#define MPI_ATTR_DELETE_OUT 34
#define MPI_ATTR_GET_IN 35
#define MPI_ATTR_GET_OUT 36
#define MPI_ATTR_PUT_IN 37
#define MPI_ATTR_PUT_OUT 38
#define MPI_COMM_COMPARE_IN 39
#define MPI_COMM_COMPARE_OUT 40
#define MPI_COMM_CREATE_IN 41
#define MPI_COMM_CREATE_OUT 42
#define MPI_COMM_DUP_IN 43
#define MPI_COMM_DUP_OUT 44
#define MPI_COMM_FREE_IN 45
#define MPI_COMM_FREE_OUT 46
#define MPI_COMM_GROUP_IN 47
#define MPI_COMM_GROUP_OUT 48
#define MPI_COMM_RANK_IN 49
#define MPI_COMM_RANK_OUT 50
#define MPI_COMM_REMOTE_GROUP_IN 51
#define MPI_COMM_REMOTE_GROUP_OUT 52
#define MPI_COMM_REMOTE_SIZE_IN 53
#define MPI_COMM_REMOTE_SIZE_OUT 54
#define MPI_COMM_SIZE_IN 55
#define MPI_COMM_SIZE_OUT 56
#define MPI_COMM_SPLIT_IN 57
#define MPI_COMM_SPLIT_OUT 58
#define MPI_COMM_TEST_INTER_IN 59
#define MPI_COMM_TEST_INTER_OUT 60
#define MPI_GROUP_COMPARE_IN 61
#define MPI_GROUP_COMPARE_OUT 62
#define MPI_GROUP_DIFFERENCE_IN 63
#define MPI_GROUP_DIFFERENCE_OUT 64
#define MPI_GROUP_EXCL_IN 65
#define MPI_GROUP_EXCL_OUT 66
#define MPI_GROUP_FREE_IN 67
#define MPI_GROUP_FREE_OUT 68
#define MPI_GROUP_INCL_IN 69
#define MPI_GROUP_INCL_OUT 70
#define MPI_GROUP_INTERSECTION_IN 71
#define MPI_GROUP_INTERSECTION_OUT 72
#define MPI_GROUP_RANK_IN 73
#define MPI_GROUP_RANK_OUT 74
#define MPI_GROUP_RANGE_EXCL_IN 75
#define MPI_GROUP_RANGE_EXCL_OUT 76
#define MPI_GROUP_RANGE_INCL_IN 77
#define MPI_GROUP_RANGE_INCL_OUT 78
#define MPI_GROUP_SIZE_IN 79
#define MPI_GROUP_SIZE_OUT 80
#define MPI_GROUP_TRANSLATE_RANKS_IN 81
#define MPI_GROUP_TRANSLATE_RANKS_OUT 82
#define MPI_GROUP_UNION_IN 83
#define MPI_GROUP_UNION_OUT 84
#define MPI_INTERCOMM_CREATE_IN 85
#define MPI_INTERCOMM_CREATE_OUT 86
#define MPI_INTERCOMM_MERGE_IN 87
#define MPI_INTERCOMM_MERGE_OUT 88
#define MPI_KEYVAL_CREATE_IN 89
#define MPI_KEYVAL_CREATE_OUT 90
#define MPI_KEYVAL_FREE_IN 91
#define MPI_KEYVAL_FREE_OUT 92
#define MPI_ABORT_IN 93
#define MPI_ABORT_OUT 94
#define MPI_ERROR_CLASS_IN 95
#define MPI_ERROR_CLASS_OUT 96
#define MPI_ERRHANDLER_CREATE_IN 97
#define MPI_ERRHANDLER_CREATE_OUT 98
#define MPI_ERRHANDLER_FREE_IN 99
#define MPI_ERRHANDLER_FREE_OUT 100
#define MPI_ERRHANDLER_GET_IN 101
#define MPI_ERRHANDLER_GET_OUT 102
#define MPI_ERROR_STRING_IN 103
#define MPI_ERROR_STRING_OUT 104
#define MPI_ERRHANDLER_SET_IN 105
#define MPI_ERRHANDLER_SET_OUT 106
#define MPI_GET_PROCESSOR_NAME_IN 107
#define MPI_GET_PROCESSOR_NAME_OUT 108
#define MPI_INITIALIZED_IN 109
#define MPI_INITIALIZED_OUT 110
#define MPI_WTICK_IN 111
#define MPI_WTICK_OUT 112
#define MPI_WTIME_IN 113
#define MPI_WTIME_OUT 114
#define MPI_ADDRESS_IN 115
#define MPI_ADDRESS_OUT 116
#define MPI_BSEND_IN 117
#define MPI_BSEND_OUT 118
#define MPI_BSEND_INIT_IN 119
#define MPI_BSEND_INIT_OUT 120
#define MPI_BUFFER_ATTACH_IN 121
#define MPI_BUFFER_ATTACH_OUT 122
#define MPI_BUFFER_DETACH_IN 123
#define MPI_BUFFER_DETACH_OUT 124
#define MPI_CANCEL_IN 125
#define MPI_CANCEL_OUT 126
#define MPI_REQUEST_FREE_IN 127
#define MPI_REQUEST_FREE_OUT 128
#define MPI_RECV_INIT_IN 129
#define MPI_RECV_INIT_OUT 130
#define MPI_SEND_INIT_IN 131
#define MPI_SEND_INIT_OUT 132
#define MPI_GET_ELEMENTS_IN 133
#define MPI_GET_ELEMENTS_OUT 134
#define MPI_GET_COUNT_IN 135
#define MPI_GET_COUNT_OUT 136
#define MPI_IBSEND_IN 137
#define MPI_IBSEND_OUT 138
#define MPI_IPROBE_IN 139
#define MPI_IPROBE_OUT 140
#define MPI_IRSEND_IN 143
#define MPI_IRSEND_OUT 144
#define MPI_ISSEND_IN 147
#define MPI_ISSEND_OUT 148
#define MPI_PACK_IN 149
#define MPI_PACK_OUT 150
#define MPI_PACK_SIZE_IN 151
#define MPI_PACK_SIZE_OUT 152
#define MPI_PROBE_IN 153
#define MPI_PROBE_OUT 154
#define MPI_RSEND_IN 157
#define MPI_RSEND_OUT 158
#define MPI_RSEND_INIT_IN 159
#define MPI_RSEND_INIT_OUT 160
#define MPI_SENDRECV_IN 163
#define MPI_SENDRECV_OUT 164
#define MPI_SENDRECV_REPLACE_IN 165
#define MPI_SENDRECV_REPLACE_OUT 166
#define MPI_SSEND_IN 167
#define MPI_SSEND_OUT 168
#define MPI_SSEND_INIT_IN 169
#define MPI_SSEND_INIT_OUT 170
#define MPI_START_IN 171
#define MPI_START_OUT 172
#define MPI_STARTALL_IN 173
#define MPI_STARTALL_OUT 174
#define MPI_TEST_IN 175
#define MPI_TEST_OUT 176
#define MPI_TESTALL_IN 177
#define MPI_TESTALL_OUT 178
#define MPI_TESTANY_IN 179
#define MPI_TESTANY_OUT 180
#define MPI_TEST_CANCELLED_IN 181
#define MPI_TEST_CANCELLED_OUT 182
#define MPI_TESTSOME_IN 183
#define MPI_TESTSOME_OUT 184
#define MPI_TYPE_COMMIT_IN 185
#define MPI_TYPE_COMMIT_OUT 186
#define MPI_TYPE_CONTIGUOUS_IN 187
#define MPI_TYPE_CONTIGUOUS_OUT 188
#define MPI_TYPE_EXTENT_IN 189
#define MPI_TYPE_EXTENT_OUT 190
#define MPI_TYPE_FREE_IN 191
#define MPI_TYPE_FREE_OUT 192
#define MPI_TYPE_HINDEXED_IN 193
#define MPI_TYPE_HINDEXED_OUT 194
#define MPI_TYPE_HVECTOR_IN 195
#define MPI_TYPE_HVECTOR_OUT 196
#define MPI_TYPE_INDEXED_IN 197
#define MPI_TYPE_INDEXED_OUT 198
#define MPI_TYPE_LB_IN 199
#define MPI_TYPE_LB_OUT 200
#define MPI_TYPE_SIZE_IN 201
#define MPI_TYPE_SIZE_OUT 202
#define MPI_TYPE_STRUCT_IN 203
#define MPI_TYPE_STRUCT_OUT 204
#define MPI_TYPE_UB_IN 205
#define MPI_TYPE_UB_OUT 206
#define MPI_TYPE_VECTOR_IN 207
#define MPI_TYPE_VECTOR_OUT 208
#define MPI_UNPACK_IN 209
#define MPI_UNPACK_OUT 210
#define MPI_WAITALL_IN 213
#define MPI_WAITALL_OUT 214
#define MPI_WAITANY_IN 215
#define MPI_WAITANY_OUT 216
#define MPI_WAITSOME_IN 217
#define MPI_WAITSOME_OUT 218
#define MPI_CART_COORDS_IN 219
#define MPI_CART_COORDS_OUT 220
#define MPI_CART_CREATE_IN 221
#define MPI_CART_CREATE_OUT 222
#define MPI_CART_GET_IN 223
#define MPI_CART_GET_OUT 224
#define MPI_CART_MAP_IN 225
#define MPI_CART_MAP_OUT 226
#define MPI_CART_SHIFT_IN 227
#define MPI_CART_SHIFT_OUT 228
#define MPI_CARTDIM_GET_IN 229
#define MPI_CARTDIM_GET_OUT 230
#define MPI_DIMS_CREATE_IN 231
#define MPI_DIMS_CREATE_OUT 232
#define MPI_GRAPH_CREATE_IN 233
#define MPI_GRAPH_CREATE_OUT 234
#define MPI_GRAPH_GET_IN 235
#define MPI_GRAPH_GET_OUT 236
#define MPI_GRAPH_MAP_IN 237
#define MPI_GRAPH_MAP_OUT 238
#define MPI_GRAPH_NEIGHBORS_IN 239
#define MPI_GRAPH_NEIGHBORS_OUT 240
#define MPI_GRAPH_NEIGHBORS_COUNT_IN 241
#define MPI_GRAPH_NEIGHBORS_COUNT_OUT 242
#define MPI_GRAPHDIMS_GET_IN 243
#define MPI_GRAPHDIMS_GET_OUT 244
#define MPI_TOPO_TEST_IN 245
#define MPI_TOPO_TEST_OUT 246
#define MPI_RECV_IDLE_IN 247
#define MPI_RECV_IDLE_OUT 248
#define MPI_CART_RANK_IN 249
#define MPI_CART_RANK_OUT 250
#define MPI_CART_SUB_IN 251
#define MPI_CART_SUB_OUT 252

//help stuff
#define MPI_COMM_MIN 39
#define MPI_COMM_MAX 60
#define MPI_GROUP_MIN 61
#define MPI_GROUP_MAX 84
#define MPI_TYPE_MIN 185
#define MPI_TYPE_MAX 208
#define MPI_ERR_MIN 95
#define MPI_ERR_MAX 106
#define MPI_CART_MIN 219
#define MPI_CART_MAX 228
#define MPI_CART2_MIN 249
#define MPI_CART2_MAX 252
#define MPI_GRAPH_MIN 233
#define MPI_GRAPH_MAX 244

//aky stuff
#define AKY_PTP_SEND 5000
#define AKY_PTP_RECV 5001

typedef struct paje_event {
  const char *name;
  const char *description;
  int id;
} s_paje_event_t, *paje_event_t;

//prototypes for aky_paje.c
void pajeDefineContainerType (const char *alias,
    const char *containerType,
    const char *name);
void pajeDefineStateType (const char *alias,
    const char *containerType,
    const char *name);
void pajeCreateContainer (double timestamp,
    const char *alias,
    const char *type,
    const char *container,
    const char *name);
void pajeDestroyContainer (double timestamp,
    const char *type,
    const char *container);
void pajeSetState (double timestamp,
    const char *container,
    const char *type,
    const char *value);
void pajePushState (double timestamp,
    const char *container,
    const char *type,
    const char *value);
void pajePopState (double timestamp,
    const char *container,
    const char *type);

void paje_header (void);
void paje_hierarchy (void);

//prototypes for aky_names.c
void name_init (void);
char *name_get (int id);

#define AKY_DEFAULT_STR_SIZE 200

#endif //__AKY_H
