#' Pair simulations and observations
#'
#' Create tables in which GETM simulations and observations are paired. The 
#' observations refer to a specific time and location, which are extracted from
#' the netcdf outputs. 
#'
#' @param ncs  character; names of the output nc file
#' @param obs character; name of the observation file (csv format)
#' @param list_vars list; list to couple the simulation variables (list names) to the names of
#'   the columns in the observation file (list values).
#' @param round_depth,round_val integer; Round depth and variable value to this many digits. No rounding if NULL.
#' @param col_date,col_x,col_y,col_depth character; column names in 'obs' file
#' @param max_time_diff Period; only select closest time if it is within this period
#'   from an observation (e.g. 'days(2)'). default = 1 day. If NULL, closest point is taken. If "exact", only
#'   exact date/time matches are used.
#' @param save_z logical; Whether to save the z-layer from which simulation value was extracted.
#'   Can be used to average over if you want to avoid multiple observations linking to the same z-layer.
#' @author
#'   Jorrit Mesman, Shajar Regev
#' @examples
#'  \dontrun{
#'    ncs = c("output_3d_20200101.nc", "output_3d_20200201.nc")
#'    obs = "obs_file.csv" # Columns "date", "x", "y", "depth", "temp", and "salinity"
#'    max_time_diff = lubridate::days(1)
#'    list_vars = list(temp = "temp",
#'                     salt = "salinity")
#'    
#'    pair_sim_and_obs(ncs = ncs, obs = obs, list_vars = list_vars)
#'  }
#' @import ncdf4 lubridate
#' @export

pair_sim_and_obs = function(ncs, obs, list_vars,
                            round_depth = NULL, round_val = NULL,
                            col_date = "date", col_x = "x", col_y = "y", col_depth = "depth",
                            save_z = TRUE,
                            max_time_diff = lubridate::days(1)){
  on.exit({
    if(exists("nc")) nc_close(nc)
  })
  
  # Load observations
  df_obs = fread(obs)
  if(any(!(unlist(list_vars) %in% names(df_obs)))){
    stop("Some of the values in 'list_vars' do not match the columns in 'obs'!")
  }
  
  if(any(df_obs[[col_depth]] > 0)){
    stop("All depths in 'df_obs' need to be negative, following definition of depth PyGETM!",
         " In the observation file it means relative to the actual water level.")
  }
  
  lst_out = vector(mode = "list", length = length(list_vars))
  names(lst_out) = names(list_vars)
  
  # First load all dates with output. This is to decide from what nc file the
  # output should be compared.
  lst_dates = lapply(ncs, function (x){
    nc = nc_open(x)
    dates = convert_nc_time(ncvar_get(nc, "time"), ncatt_get(nc, "time"))
    nc_close(nc)
    data.table(date = dates,
               nc = x)
  })
  df_dates = rbindlist(lst_dates)
  setorder(df_dates, date)
  if(any(duplicated(df_dates$date))) stop("Multiple nc files with same date!")
  
  # Add nc file with closest observation, incl. max_time_diff
  df_dates[, sim_date := date]
  df_obs = df_dates[df_obs, roll = "nearest", on = "date"]
  
  if(identical(max_time_diff, "exact")){
    df_obs[!(date == sim_date), nc := NA]
  }else if(!is.null(max_time_diff)){
    df_obs[abs(difftime(date, sim_date)) > max_time_diff, nc := NA]
  }
  
  ncs_to_read = unique(df_obs[!is.na(nc), nc])
  progress = 0
  pb = txtProgressBar(min = progress, max = length(ncs_to_read) * length(list_vars), style = 3)
  
  for(i in seq_along(ncs_to_read)){
    # Open NC
    nc = nc_open(ncs_to_read[i])
    if(any(!(names(list_vars) %in% names(nc$var)))){
      stop("Some of the names in 'list_vars' do not match the variable names in ",
           ncs_to_read[i] , "!")
    }
    
    nc_dates = convert_nc_time(ncvar_get(nc, "time"), ncatt_get(nc, "time"))
    
    # Cut obs file to relevant dates
    df_obs_i = df_obs[nc == ncs_to_read[i]]
    
    # Load x and y dimensions
    x_dim = ncvar_get(nc, "xt")[, 1] # Equidistant grid, so column 1 is the same as any other
    y_dim = ncvar_get(nc, "yt")[1,]
    
    # Load heights of cell interfaces
    m_zft = ncvar_get(nc, varid = "zft")
    if(length(dim(m_zft)) == 3L) dim(m_zft) = c(dim(m_zft), 1)
    
    for(j in names(list_vars)){
      if(all(is.na(df_obs_i[, get(list_vars[[j]])]))){
        progress = progress + 1
        setTxtProgressBar(pb, progress)
        next
      }
      
      m_var = ncvar_get(nc, varid = j)
      if(length(dim(m_var)) == 3L) dim(m_var) = c(dim(m_var), 1)
      
      df_sim = slice_matrix_with_obs(mtrx = m_var, mtrx_zft = m_zft,
                                     x_dim = x_dim, y_dim = y_dim,
                                     time_dim = nc_dates,
                                     observations = df_obs_i[!is.na(get(list_vars[[j]])),
                                                             c(col_date, "sim_date", col_x, col_y, col_depth, list_vars[[j]]), with = F],
                                     col_depth = col_depth, sim_date = "sim_date",
                                     col_x = col_x, col_y = col_y,
                                     save_z = save_z)
      # Store info in list
      lst_out[[j]][[length(lst_out[[j]]) + 1]] = df_sim
      
      progress = progress + 1
      setTxtProgressBar(pb, progress)
    }
    
    nc_close(nc)
    rm(nc) # Has to be done to play nice with the on.exit statement
  }
  
  return(lapply(lst_out, rbindlist))
}
