/**
 * Copyright 2018
 *
 * \author   <scott.wales@unimelb.edu.au>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * BLAH
 */
#include "blah.h"
#include "netcdf.h"
#include <argp.h>
#include <assert.h>
#include <glob.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "async.h"
#include "error.h"
#include "read_chunked.h"

#define TAG_DECOMP 1
#define TAG_CHUNK 2

struct args_t {
  const char *output;
  int deflate_level;
  int shuffle;
  bool force;
  bool remove;
};

// Copy a (possibly collated) field in NetCDF mode
void copy_netcdf(int ncid_out, int varid_out, int ncid_in, int varid_in) {
  int ndims;
  NCERR(nc_inq_varndims(ncid_in, varid_in, &ndims));

  size_t in_offset[ndims];
  size_t out_offset[ndims];
  size_t local_size[ndims];
  size_t total_size[ndims];

  get_collation_info(ncid_in, varid_in, out_offset, local_size, total_size,
                     ndims);

  size_t size = 1;
  for (int d = 0; d < ndims; ++d) {
    size *= local_size[d];
    in_offset[d] = 0;
  }
  if (size > (size_t)2<<30) {
    CERR(-1, "Variable size too large");
  }

  // Enough size for float64_t
  void *buffer = malloc(size * 8);
  NCERR(nc_get_vara(ncid_in, varid_in, in_offset, local_size, buffer));
  NCERR(nc_put_vara(ncid_out, varid_out, out_offset, local_size, buffer));
  free(buffer);
}

// Copy NetCDF attributes, either globally or for a variable
void copy_attrs(int ncid_out, int varid_out, int ncid_in, int varid_in,
                int natts) {
  size_t buffer_len = 1024 * 8;
  char buffer[buffer_len];

  for (int a = 0; a < natts; ++a) {
    char attname[NC_MAX_NAME + 1];
    nc_type atttype;
    size_t attlen;
    NCERR(nc_inq_attname(ncid_in, varid_in, a, attname));
    NCERR(nc_inq_att(ncid_in, varid_in, attname, &atttype, &attlen));
    if (varid_in == NC_GLOBAL && strcmp(attname, "NumFilesInSet") == 0) {
      continue; // Don't copy
    }
    // Check the buffer is big enough
    assert(attlen * 8 < buffer_len);
    NCERR(nc_get_att(ncid_in, varid_in, attname, buffer));
    if (varid_in == NC_GLOBAL && strcmp(attname, "filename") == 0) {
      // Replace filename attribute with new collated filename
      NCERR(nc_inq_path(ncid_out, &attlen, NULL));
      // Check the buffer is big enough
      assert(attlen * 8 < buffer_len);
      NCERR(nc_inq_path(ncid_out, NULL, buffer));
    }
    NCERR(nc_put_att(ncid_out, varid_out, attname, atttype, attlen, buffer));
  }
}

// Copy NetCDF headers and uncollated variables from file at in_path to file at
// out_path
void init(const char *in_path, const char *out_path,
          const struct args_t *args) {
  int in_file;
  int out_file;
  int err;

  // Open both files
  NCERR(nc_open(in_path, NC_NOWRITE, &in_file));

  int out_flags = NC_NETCDF4;
  if (!args->force) {
    out_flags |= NC_NOCLOBBER;
  }
  err = nc_create(out_path, out_flags, &out_file);
  if (err == NC_EEXIST) {
    log_message(LOG_ERROR,
                "ERROR: Output file '%s' already exists (try --force)",
                out_path);
    MPI_Abort(MPI_COMM_WORLD, err);
  } else {
    NCERR(err);
  }

  int ndims;
  int nvars;
  int natts;
  NCERR(nc_inq(in_file, &ndims, &nvars, &natts, NULL));

  // Copy global attributes
  copy_attrs(out_file, NC_GLOBAL, in_file, NC_GLOBAL, natts);

  // Get count of unlimited dimensions
  err = nc_inq_unlimdims(in_file, NULL, NULL);
  int nunlimdim;
  if (err != NC_ENOTNC4) {
    NCERR(nc_inq_unlimdims(in_file, &nunlimdim, NULL));
  } else {
    nunlimdim = 1;
  }

  // Get the dimension ids
  int unlimdims[nunlimdim];
  if (err != NC_ENOTNC4) {
    NCERR(nc_inq_unlimdims(in_file, NULL, unlimdims));
  } else {
    NCERR(nc_inq_unlimdim(in_file, unlimdims));
    if (unlimdims[0] < 0) {
      nunlimdim = 0;
    }
  }

  // Copy dimensions
  for (int d = 0; d < ndims; ++d) {
    char name[NC_MAX_NAME + 1];
    size_t len;
    int dimid;
    int varid;

    NCERR(nc_inq_dim(in_file, d, name, &len));
    log_message(LOG_DEBUG, "checking variable %s", name);

    // Check if the variable with the same name is collated
    err = nc_inq_varid(in_file, name, &varid);
    if (err == NC_ENOTVAR) {
      // No associated variable
      continue;
    }
    NCERR(err);

    if (is_collated(in_file, varid)) {
      // If so get the full length
      len = get_collated_dim_len(in_file, name);
    }

    for (int ud = 0; ud < nunlimdim; ++ud) {
      if (d == unlimdims[ud]) {
        len = NC_UNLIMITED;
      }
    }

    // Create the out dim
    NCERR(nc_def_dim(out_file, name, len, &dimid));
    assert(dimid == d);
  }

  // Copy variables
  for (int v = 0; v < nvars; ++v) {
    char name[NC_MAX_NAME + 1];
    nc_type type;
    int ndims;
    int natts;

    NCERR(nc_inq_var(in_file, v, name, &type, &ndims, NULL, &natts));

    int dimids[ndims];
    NCERR(nc_inq_vardimid(in_file, v, dimids));

    log_message(LOG_DEBUG, "Defining variable %s", name);

    int out_v;
    NCERR(nc_def_var(out_file, name, type, ndims, dimids, &out_v));

    // Set up chunking
    int storage;
    size_t chunk[ndims];
    NCERR(nc_inq_var_chunking(in_file, v, &storage, chunk));
    if (storage == NC_CHUNKED) {
      NCERR(nc_def_var_chunking(out_file, out_v, storage, chunk));
    }

    // Set up compression
    int shuffle;
    int deflate;
    int deflate_level;
    err = nc_inq_var_deflate(in_file, v, &shuffle, &deflate, &deflate_level);
    if (err == NC_ENOTNC4) {
        shuffle = 1;
        deflate = 1;
        deflate_level = 6;
    } else {
        NCERR(err);
    }

    // Option to override compression
    if (args->deflate_level != -1) {
      deflate_level = args->deflate_level;
      deflate = 1;
    }
    if (args->shuffle != -1) {
      shuffle = args->shuffle;
    }

    if (shuffle || deflate) {
      NCERR(
          nc_def_var_deflate(out_file, out_v, shuffle, deflate, deflate_level));
    }

    // Copy attributes
    copy_attrs(out_file, v, in_file, v, natts);

    // If the field is not collated copy it now
    if (!is_collated(in_file, v)) {
      log_message(LOG_INFO, "Uncollated NetCDF copy of %s", name);
      copy_netcdf(out_file, out_v, in_file, v);
    }
  }

  // Close the files
  NCERR(nc_close(in_file));
  NCERR(nc_close(out_file));
}

// Copy contiguous variables - no chunking means no compression, so we can just
// use NetCDF
void copy_contiguous(const char *out_path, char **in_paths, int n_in) {

  int out_nc4;
  NCERR(nc_open(out_path, NC_WRITE, &out_nc4));

  int nvars;
  NCERR(nc_inq_nvars(out_nc4, &nvars));

  for (int i = 0; i < n_in; ++i) {
    int in_nc4;
    NCERR(nc_open(in_paths[i], NC_NOWRITE, &in_nc4));

    for (int v = 0; v < nvars; ++v) {
      // Check the metadata of the output file to see if this variable is
      // contiguous
      char varname[NC_MAX_NAME + 1];
      NCERR(nc_inq_varname(out_nc4, v, varname));
      int storage;
      NCERR(nc_inq_var_chunking(in_nc4, v, &storage, NULL));

      if (storage == NC_CONTIGUOUS && is_collated(out_nc4, v)) {
        log_message(LOG_INFO, "NetCDF copy of %s from %s", varname,
                    in_paths[i]);
        copy_netcdf(out_nc4, v, in_nc4, v);
      }
    }

    NCERR(nc_close(in_nc4));
  }

  NCERR(nc_close(out_nc4));
}

void file_match_check(bool test, const char *filea, const char *fileb,
                      const char *message) {
  if (!test) {
    fprintf(stderr, "ERROR: %s <%s> <%s>\n", message, filea, fileb);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
}

void add_metadata(const char *file, int argc, char **const argv) {
  int ncid;
  NCERR(nc_open(file, NC_WRITE, &ncid));

  // Remove decomposition info from collated variables
  int nvars;
  NCERR(nc_inq_nvars(ncid, &nvars));
  for (int v = 0; v < nvars; ++v) {
    int err = nc_del_att(ncid, v, "domain_decomposition");
    if (err == NC_ENOTATT) {
      continue;
    }
    NCERR(err);
  }

  // Build up a string with the command line
  size_t argv_len[argc];
  size_t total_len = 0;
  for (int a = 0; a < argc; ++a) {
    argv_len[a] = strlen(argv[a]) + 1;
    total_len += argv_len[a];
  }
  char *commandline = malloc(total_len);
  size_t offset = 0;
  for (int a = 0; a < argc; ++a) {
    memcpy(commandline + offset, argv[a], argv_len[a]);
    offset += argv_len[a];
    commandline[offset - 1] = ' ';
  }
  commandline[total_len - 1] = '\0';

  if (total_len > 100) {
    log_message(LOG_WARNING,
                "WARNING: Output file history is long, "
                "consider escaping wildcards like `./mppnccombine-fast "
                "input.nc.\\* -o output.nc`");
  }

  // Get the current time
  time_t now = time(NULL);
  struct tm *time_tm = gmtime(&now);

  char timestamp[4 + 2 + 2 + 1 + 2 + 2 + 3];
  int err = strftime(timestamp, sizeof(timestamp), "%Y%m%dT%H%MZ:", time_tm);
  if (err == 0) {
    CERR(-1, "Error calculating timestamp");
  }

  // Get the old history
  size_t old_history_len;
  err = nc_inq_attlen(ncid, NC_GLOBAL, "history", &old_history_len);
  if (err == NC_ENOTATT) {
    log_message(LOG_DEBUG, "No previous history");
    old_history_len = 0;
  }

  size_t history_len = old_history_len + sizeof(timestamp) + total_len;
  log_message(LOG_DEBUG, "History length %zu", history_len);
  char *history = malloc(history_len);

  if (old_history_len > 0) {
    NCERR(nc_get_att_text(ncid, NC_GLOBAL, "history", history));
    log_message(LOG_DEBUG, "%s", history);
    history[old_history_len - 1] = '\n';
  }

  memcpy(history + old_history_len, timestamp, sizeof(timestamp));
  log_message(LOG_DEBUG, "%s", history);
  history[old_history_len + sizeof(timestamp) - 1] = ' ';
  memcpy(history + old_history_len + sizeof(timestamp), commandline, total_len);
  log_message(LOG_DEBUG, "%s", history);

  NCERR(nc_put_att_text(ncid, NC_GLOBAL, "history", history_len, history));

  free(history);
  free(commandline);
  NCERR(nc_close(ncid));
}

static char doc[] =
    "\nQuickly collate MOM model files\n\nGathers the INPUT MOM model files "
    "(provided e.g. with a shell glob) and joins them along their horizontal "
    "dimensions into a single NetCDF file";

static struct argp_option opts[] = {
    {"INPUT", 0, 0, OPTION_DOC, "Input NetCDF file", 0},
    {"output", 'o', "OUTPUT", OPTION_NO_USAGE, "Output NetCDF file", 1},
    {"force", 'f', 0, 0, "Combine even if output file present", 2},
    {"remove", 'r', 0, 0, "Remove the input files after completion", 3},
    {"deflate", 'd', "[0-9]", 0, "Override compression level (slower)", 4},
    {"shuffle", -1, 0, 0, "Force enable shuffle filter (slower)", 5},
    {"no-shuffle", -2, 0, 0, "Force disable shuffle filter (slower)", 6},
    {"verbose", 'v', 0, 0, "Be verbose", 7},
    {"quiet", 'q', 0, 0, "Be quiet (no warnings)", 8},
    {"debug", -3, 0, 0, "Debug info", 9},
    {"help", '?', 0, 0, "Print this help list", 10},
    {0},
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  struct args_t *args = state->input;
  int rank;
  char *endptr;

  switch (key) {
  case 'o':
    args->output = arg;
    break;
  case 'd':
    args->deflate_level = strtol(arg, &endptr, 10);
    if (*endptr != '\0') {
      CERR(-1, "Bad deflate value");
    }
    if (args->deflate_level < 0) {
      CERR(-1, "Bad deflate value");
    }
    if (args->deflate_level > 9) {
      CERR(-1, "Bad deflate value");
    }
    break;
  case -1:
    args->shuffle = 1;
    break;
  case -2:
    args->shuffle = 0;
    break;
  case 'f':
    args->force = true;
    break;
  case 'r':
    args->remove = true;
    break;
  case 'v':
    set_log_level(LOG_INFO);
    break;
  case 'q':
    set_log_level(LOG_ERROR);
    break;
  case -3:
    set_log_level(LOG_DEBUG);
    break;
  case '?':
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      argp_state_help(state, stderr,
                      ARGP_HELP_SHORT_USAGE | ARGP_HELP_DOC | ARGP_HELP_LONG);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, 0);
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }

  return 0;
}

static struct argp argp = {
    .parser = parse_opt,
    .options = opts,
    .doc = doc,
    .args_doc = "--output=OUTPUT INPUT [INPUT ...]",
};

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  set_log_level(LOG_WARNING);

  int comm_rank;
  int comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int current_file_idx = 0;
  MPI_Win current_file_win;
  MPI_Win_create(&current_file_idx, sizeof(current_file_idx),
                 sizeof(current_file_idx), MPI_INFO_NULL, MPI_COMM_WORLD,
                 &current_file_win);

  int arg_index;
  struct args_t args = {0};
  args.deflate_level = -1;
  args.shuffle = -1;
  args.force = false;
  args.remove = false;

  unsigned int argp_flags = ARGP_NO_EXIT | ARGP_NO_HELP;
  if (comm_rank != 0) {
    argp_flags |= ARGP_SILENT;
  }

  error_t argp_error =
      argp_parse(&argp, argc, argv, argp_flags, &arg_index, &args);
  if (comm_rank == 0) {
    if (args.output == NULL) {
      log_message(LOG_ERROR, "No output file specified");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (arg_index == argc) {
      log_message(LOG_ERROR, "No input files specified");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (comm_size < 2) {
      log_message(LOG_ERROR, "Please run with at least 2 MPI processes");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (argp_error != 0) {
      log_message(LOG_ERROR, "ERROR parsing arguments");
      MPI_Abort(MPI_COMM_WORLD, argp_error);
    }
  }

  const char *out_path = args.output;

  // Glob inputs
  int err = -1;
  glob_t globs = {0};
  int glob_flags = GLOB_BRACE | GLOB_TILDE;
  for (int i = arg_index; i < argc; ++i) {
    err = glob(argv[i], glob_flags, NULL, &globs);
    glob_flags |= GLOB_APPEND;
  }

  if ((err != 0 || globs.gl_pathc < 1) && comm_rank == 0) {
    log_message(LOG_ERROR, "No matching input files found");
    MPI_Abort(MPI_COMM_WORLD, -1);
    exit(-1); // Tell clang-tidy we're bailing on glob error
  }

  int writer_rank = 0;

  if (comm_rank == writer_rank) {
    // Copy metadata and un-collated variables
    fprintf(stdout, "\nCopying non-collated variables\n");
    init(globs.gl_pathv[0], out_path, &args);
    // Copy contiguous variables using NetCDF
    fprintf(stdout, "\nCopying contiguous variables\n");
    copy_contiguous(out_path, globs.gl_pathv, globs.gl_pathc);

    fprintf(stdout, "\nCopying chunked variables\n");
    double t_start = MPI_Wtime();
    size_t total_size = run_async_writer(out_path);
    double t_end = MPI_Wtime();
    double total_size_gb = total_size / pow(1024, 3);

    // Clean up attributes & add to history
    add_metadata(out_path, argc, argv);

    fprintf(stdout,
            "\nTotal compressed size %.2f GiB | Time %.2fs | %.2f MiB / sec\n",
            total_size_gb, t_end - t_start,
            total_size / pow(1024, 2) / (t_end - t_start));
  } else {
    int increment = 1;
    unsigned int my_file_idx = -1;

    log_message(LOG_DEBUG, "Starting read");

    // Atomic post-addition of increment to current_file_idx
    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, current_file_win);
    MPI_Fetch_and_op(&increment, &my_file_idx, MPI_INT, 0, 0, MPI_SUM,
                     current_file_win);
    MPI_Win_unlock(0, current_file_win);

    while (my_file_idx < globs.gl_pathc) {
      // Read chunked variables using HDF5, sending data to the async_writer to
      // be written
      log_message(LOG_DEBUG, "%d %s", my_file_idx, globs.gl_pathv[my_file_idx]);
      copy_chunked(globs.gl_pathv[my_file_idx], writer_rank);

      MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, current_file_win);
      MPI_Fetch_and_op(&increment, &my_file_idx, MPI_INT, 0, 0, MPI_SUM,
                       current_file_win);
      MPI_Win_unlock(0, current_file_win);
    }
    close_async(0);
    log_message(LOG_DEBUG, "Finished read");
  }

  if (comm_rank == writer_rank && args.remove) {
    log_message(LOG_INFO, "Cleaning inputs");
    for (unsigned int my_file_idx = 0; my_file_idx < globs.gl_pathc; ++my_file_idx) {
      log_message(LOG_DEBUG, "unlink %s", globs.gl_pathv[my_file_idx]);
      unlink(globs.gl_pathv[my_file_idx]);
    }
  }

  log_message(LOG_DEBUG, "Cleanup glob");
  globfree(&globs);
  log_message(LOG_DEBUG, "Cleanup win");
  // MPI_Win_free(&current_file_win);

  log_message(LOG_DEBUG, "MPI_Finalize");
  return MPI_Finalize();
}
