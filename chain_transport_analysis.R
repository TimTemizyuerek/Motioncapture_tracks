###################################################################################################################################################
##
## Project: Chain transport
## Script purpose: Preliminary motion capture analysis
## Start date: September 2022
## Author: Tim Temizyürek (timtemiz@gmail.com)
##
###################################################################################################################################################

## file patch, libraries & custom functions ####

     ## paths
     dir_data = "C:/Users/timte/Desktop/Konstanz/Chain transport/pfolder/"
     # dir_github = "C:/Users/timte/Desktop/Konstanz/Chain transport/Github repository/"
     
     ## custom functions
     
     ## make raw dataframe into concise matrix (I am speed)
     ## produces: matrix with cols: frame_number, time_in_deciseconds, X, Y, Z & ID
     m.concise.dataframe = function(input_ma) {
          
          ## prepare output matrix
          concise_matrix = matrix(nrow=0, ncol=6)
          
          ## create ID counter
          ID_counter = 0
          
          ## reorganize dataframe
          for (n in seq(from=3, to=(ncol(input_ma)-2), by=3)){
               
               ## add to ID counter
               ID_counter = ID_counter + 1
               
               ## extract coordinates for one marker
               runner_m = as.matrix(input_ma[,c(1:2,n:(n+2))])
               
               ## find and extract rows without NAs 
               good_subset = runner_m[which(!is.na(runner_m[,3])),]
               
               ## generate matrix
               runner_merger = matrix(c(good_subset[,1], 
                                        good_subset[,2], 
                                        good_subset[,3],
                                        good_subset[,4],
                                        good_subset[,5],
                                        rep(ID_counter,nrow(good_subset))), ncol = 6,
                                        )
               
               ## bind dataframe
               concise_matrix = rbind(concise_matrix, runner_merger)
          
               ## progress update
               print(paste("concise: ",round(n/(ncol(input_ma)),5)*100," %",sep=""))
          }
          
          ## name the cols of the matrix for convenience
          colnames(concise_matrix) = c("frame_number", 
                                       "time_in_deciseconds",
                                       "X",
                                       "Y",
                                       "Z",
                                       "ID")
          
          ## return the concise matrix
          return(concise_matrix)
          
     }
          
     ## takes pretty data and extracts the first and last of each tracklet in two separate matrices
     ## produces: first- and last-data matrix from a pretty data matrix
     first.last.finder = function (pretty_data) {
          
          ## store rownumbers
          first_index = rep(NA, length(unique(pretty_data[,"ID"])))
          last_index = rep(NA, length(unique(pretty_data[,"ID"])))          
          
          ## identify first and last datapoint of each tracklet
          for (n in 1:length(unique(pretty_data[,"ID"]))) {
               
               ## run through IDs
               runner_ID = unique(pretty_data[,"ID"])[n]
               
               ## take last and first marker for each ID
               first_index[n] = head(which(pretty_data[,"ID"] == runner_ID),1)
               last_index[n] = tail(which(pretty_data[,"ID"] == runner_ID),1)
               
               ## progress update
               print(paste(round(n/(length(unique(pretty_data[,"ID"]))),5)*100," %",sep=""))
          }
          
          ## extract first and last datapoints
          stitching_data_first = as.matrix(pretty_data[first_index,, drop=FALSE])
          stitching_data_last = as.matrix(pretty_data[last_index,, drop=FALSE])
          
          
          ## create return list
          resturn_list = list(stitching_data_first, stitching_data_last)
          
          ## output
          return(resturn_list)
     }
     
     ## finds candidates for fragment stitching
     ## produces: first- and last-data matrix, list with stitched tracks
     fragment.stitcher = function (u_first, u_last, u_track_ID, u_time_window, u_x_diff, u_y_diff) {
       
          ## for testing
          # u_first = seq_results[[1]]; u_last = seq_results[[2]]; u_track_ID = seq_results[[3]]; u_time_window = 10; u_x_diff = 10; u_y_diff = 10
          
          ## variables for the while loop
          proceed = TRUE
          n=1
          
          ## loop over the endpoints
          while (proceed == TRUE) {
               
               ## take endpoint
               runner_endpoint = u_last[n,,drop=FALSE]

               ## select all startpoints, that the come after runner_endpoint and are within u_time_window
               runner_candidates_time = u_first[which(u_first[,"frame_number"] > as.numeric(runner_endpoint[,"frame_number"]) & 
                                                      u_first[,"frame_number"] <= (as.numeric(runner_endpoint[,"frame_number"])+u_time_window)),, drop=FALSE]
               
               ## select all rows differ less than u_x_diff units on X
               runner_candidate_time_x = runner_candidates_time[which((abs(runner_candidates_time[,"X"] - as.numeric(runner_endpoint[,"X"]))) <= u_x_diff),, drop=FALSE]
               
               ## select all rows differ less than u_y_diff units on Y
               runner_candidate_time_x_y = runner_candidate_time_x[which((abs(runner_candidate_time_x[,"Y"] - as.numeric(runner_endpoint[,"Y"]))) <= u_y_diff),, drop=FALSE]
               
               ## if one tracklet left: stitch
               if (nrow(runner_candidate_time_x_y) == 1) {
                    
                    ## update u_track_ID
                    runner_IDs = as.numeric(c(u_track_ID[[runner_endpoint[,"ID"]]], runner_endpoint[,"ID"], runner_candidate_time_x_y[,"ID"]))
                    u_track_ID[[min(runner_IDs)]] = as.numeric(c(u_track_ID[[runner_endpoint[,"ID"]]], runner_endpoint[,"ID"], runner_candidate_time_x_y[,"ID"]))
                         
                    ## remove endpoint of runner_endpoint
                    u_last = u_last[-which(u_last[,"ID"] == runner_endpoint[,"ID"]),, drop=FALSE]
                    
                    ## update new endpoint ID
                    u_last[which(u_last[,"ID"] == runner_candidate_time_x_y[,"ID"]),"ID"] = as.numeric(runner_endpoint[,"ID"])
                    
                    ## remove now stitched startpoint runner_candidate_time_x_y
                    u_first = u_first[-which(u_first[,"ID"] == runner_candidate_time_x_y[,"ID"]),, drop=FALSE]
                    
               } else {
                    
                    ## update u_track_ID
                    # u_track_ID[[runner_endpoint[,"ID"]]] = as.numeric(c(u_track_ID[[runner_endpoint[,"ID"]]], runner_endpoint[,"ID"]))
                    
                    ## go to next tracklet
                    n = n + 1
                    
                    ## stop if there is no more next tracklet
                    proceed = ifelse(n > length(u_last[,"ID"]), FALSE, TRUE)
               }
          
          }
          
          ## only keep unique items in u_track_ID
          u_track_ID = lapply(u_track_ID, unique)
          
          ## make return list
          return_list = list(u_first, u_last, u_track_ID)     
          
          ## function output
          return(return_list)
     }
          
     ## sequentially apply fragment.stitcher to the data, updating parameters on the way
     ## produces: first- and last-data matrix, list with stitched tracks
     sequential.fragment.stitcher = function(u_first, u_last, u_track_ID, u_time_window_start, u_time_window_steps,u_x_diff_start, u_x_diff_steps,u_y_diff_start, u_y_diff_steps,iterations) {
          
          ## run the function on the starting values
          runner_data = fragment.stitcher(u_first, u_last, u_track_ID, u_time_window_start, u_x_diff_start, u_y_diff_start)
          
          if (iterations >= 1) {
               for (n in 1:iterations){
                    
                    ## print progress
                    print(paste(round((n/iterations)*100,3)," %", sep=""))
                    
                    ## update parameters
                    time_window = u_time_window_start + (n*u_time_window_steps) 
                    x_diff = u_x_diff_start + (n*u_x_diff_steps) 
                    y_diff = u_y_diff_start + (n*u_y_diff_steps)
                    
                    ## rerun function
                    runner_data = fragment.stitcher(runner_data[[1]], runner_data[[2]], runner_data[[3]], time_window, x_diff, y_diff)
               }
          } else if (iterations < 1) {
               
               stop("ERROR: iterations must be 1 or higher")
               
          } 
          
          return(runner_data)
     }
     
     ## extract pretty data from stitching
     ## produces a list, where each entry is a stitched track
     track.extractor = function (u_track_ID) {

          ## list of matches
          matches_list = vector(mode='list', length=length(u_track_ID))
               
          ## run through loop
          for (n in 1:length(u_track_ID)){
               
               ## go through track IDs
               while_runner = u_track_ID[[n]]
               
               ## use index to find matching elements in list
               proceed = TRUE
               while(proceed == TRUE){
                    
                    ## take length for below
                    pre_length = length(while_runner)
                    
                    ## use index to find other elements
                    while_runner = unique(unlist(u_track_ID[while_runner]))
                    
                    ## end loop if nothing new was added in last iteration
                    if (pre_length == length(while_runner)) {proceed = FALSE}
               }
               
               ## delete elements from list
               u_track_ID[sort(while_runner)[2:length(while_runner)]] = NA
               # u_track_ID[(while_runner)] = NA
               
               ## add stitched tracklets to output list
               if (is.null(while_runner)) {matches_list[[n]] = while_runner} else {matches_list[[min(while_runner)]] = while_runner}
          }
          
          ## format output
          output_list = matches_list[!sapply(matches_list, is.null)]
          return(output_list)
     }
     
     ## assembles raw data tracks from track.extractor
     ## produces: concise df with stitched tracks
     track.data = function (u_track_ID, u_data) {
          
          ## add ID_track to df
          u_data = cbind(u_data, rep(NA, nrow(u_data)))
          colnames(u_data)[ncol(u_data)] = "ID_track"
          
          ## output list for dataframes containing track data
          track_list = vector(mode='list', length=length(u_track_ID))
          
          ## run through tracks
          for (n in 1:length(u_track_ID)) {
               
               ## run through list
               runner = u_track_ID[[n]]
               
               ## assemble df, add ID_track col put into list
               runner_df = u_data[which(u_data[,"ID"] %in% runner),]
               runner_df[,"ID_track"] = n
               track_list[[n]] = runner_df
          }
          
          ## output
          return(track_list)
     }
     
     ## computes window based on parameters
     ## produces: dataframe with all window sizes
     windows_from_parameter = function(time_start, time_steps, x_start, x_steps, y_start, y_steps, iterations) {
          
          ## create dataframe to store windows
          window_df = data.frame(matrix(ncol = 4, nrow = iterations))
          colnames(window_df) <- c("time", "x","y","iteration")
          
          ## calculate window in all iterations
          for (n in 1:iterations) {
               
               ## calculate window
               window_df[n,"time"] = time_start + (time_steps*(n-1))
               window_df[n,"x"] = x_start + (x_steps*(n-1))
               window_df[n,"y"] = y_start + (y_steps*(n-1))
               window_df[n,"iteration"] = n
          }
          
          return(window_df)
          
     }
     
     ## computes distances between datapoins, within tracklets
     ## produces: adds a distance column to u_data
     distance.calculator = function (u_data) {
          
          ## add distance col
          u_data = cbind(u_data, rep(NA, nrow(u_data)))
          colnames(u_data)[ncol(u_data)] = "distance"
                                     
          ## calculate all distances
          for (n in 1:(nrow(u_data)-1)) {
          
               ## calculate distance
               u_data[n,"distance"] = round( sqrt((as.numeric(u_data[n,"X"])-as.numeric(u_data[n+1,"X"]))^2+(as.numeric(u_data[n,"Y"])-as.numeric(u_data[n+1,"Y"]))^2),3)
          
               ## timer
               print(paste("distance: ", round(n/nrow(u_data),3), " %", sep=""))
          }
          
          ## remove distances calculated between tracklets
          u_data[which(diff(u_data[,"ID"]) != 0),"distance"] = NA
          
          ## return dataframe
          return(u_data)
          
     }
     
     ## takes a track and identifies waiting times
     ## produce: list matching track_data with waiting times
     wait.calculator = function(u_track_data, u_min_waiting_time) {
          
          ## waiting times list
          waiting_times = vector(mode='list', length=length(u_track_data))
          
          for (n in 1:length(u_track_data)){
               
               ## run through tracks
               runner_track = u_track_data[[n]]
          
               ## find places where distances don't go above 1 (which is measurement error)
               runner_stationary = runner_track[which(runner_track[,"distance"] <= 1.5),,drop=FALSE]
               
               ## identify breaks in stationary phases
               runner_breaks = which(diff(runner_stationary[,"frame_number"]) == 1)
               
               ## split potential waiting times & remove too short ones
               split_candidates = split(runner_breaks, cumsum(c(1, diff(runner_breaks) != 1)))
               split_candidates = split_candidates[which(lapply(split_candidates, length) >= u_min_waiting_time)]
               
               ## re-extract data
               candidates_data = lapply(split_candidates, function(x) runner_track[x,,drop=FALSE])
               
               plot(candidates_data[[1]][,"X"], candidates_data[[1]][,"Y"])
               
               
               
               
               
               ## find lengths of waiting times
               runner_waiting_times = diff(runner_breaks)
               
               ## remove short waiting times and put in return list
               waiting_times[[n]] = runner_waiting_times[which(runner_waiting_times >= u_min_waiting_time)]
               
          }
           
          ## output
          return(waiting_times)
          
     }
          
     
## 1. DATA PREPARATION (LOADING, TRANSFORMATION, STITCHING) ####
     ## load raw data and transform to concise matrix ####
     
     ## raw_data = read.table(file = paste(dir_data, "Quality_check0001_04102022.tsv", sep=""), sep = '\t', header=FALSE, skip=12, fill=FALSE)
     ## raw_data = read.table(file = paste(dir_data, "07102022_data_collection_0001.tsv", sep=""), sep = '\t', header=FALSE, skip=12, fill=FALSE)
     
     ## automate for several files
     tsv_files = list.files(dir_data, pattern=".tsv")
     
     ## loop through files
     for (n in 1:length(tsv_files)) {
          
          ## print process
          print(paste("processing file ", n, " out of", length(tsv_files)), sep="")
          
          ## load file
          runner_raw_data = read.table(file = paste(dir_data, tsv_files[n], sep=""), sep = '\t', header=FALSE, skip=12, fill=FALSE)
          
          ## remove "measured-col" from the dataframe (if necessary) and make matrix
          col_to_remove = as.numeric(which(apply(runner_raw_data, 2, function(x) sum(unique(x) == "Measured", na.rm=TRUE)) == 1))
          if (sum(col_to_remove) != 0) {runner_raw_data = as.matrix(runner_raw_data[,-col_to_remove])}
          
          ## make concise data
          runner_raw_concise_data = m.concise.dataframe(runner_raw_data)
          
          ## remove objects to free up space
          rm(runner_raw_data)
          
          ## calculate distances 
          runner_raw_concise_data = distance.calculator(runner_raw_concise_data)
          
          ## remove tracklets with too high/low Z values
          ## use Z values of the trail outline + 10 cm buffer
          raw_trail_outline_data = read.table(file = paste(dir_data, "Trail_outline.tsv", sep=""), sep = '\t', header=FALSE, skip=12, fill=FALSE, nrows=2)
          clean_trail_data = m.concise.dataframe(raw_trail_outline_data); clean_trail_data = clean_trail_data[which(clean_trail_data[,"frame_number"] == 1),]; clean_trail_data = clean_trail_data[-78,]
          ## generate cutoff values
          low_cutoff_z = min(clean_trail_data[,"Z"]) - 100; high_cutoff_z = max(clean_trail_data[,"Z"]) + 100
          ## extract affected IDs
          IDs_to_remove_because_z = unique(runner_raw_concise_data[which(runner_raw_concise_data[,"Z"] <= low_cutoff_z | runner_raw_concise_data[,"Z"] >= high_cutoff_z),][,"ID"])
          ## remove affected IDs
          if (length(IDs_to_remove_because_z) > 0) {runner_raw_concise_data = runner_raw_concise_data[-which(runner_raw_concise_data[,"ID"] %in% IDs_to_remove_because_z),]}
          
          ## remove tracklets with big jumps in Z values
          ## calculate differences
          z_distances = diff(runner_raw_concise_data[,"Z"])
          ## remove differences between different tracklets
          z_distances[which(diff(runner_raw_concise_data[,"ID"]) != 0)] = NA
          ## identify tracklets with jumps of more than 5cm
          IDs_to_remove_because_of_z_diff = runner_raw_concise_data[which(z_distances >= 50),"ID"]
          ## remove them
          if (length(IDs_to_remove_because_of_z_diff) > 0) {runner_raw_concise_data = runner_raw_concise_data[-which(runner_raw_concise_data[,"ID"] %in% IDs_to_remove_because_of_z_diff),]}
          
          ## export new dataframe to txt
          write.table(runner_raw_concise_data, file=paste(dir_data, paste("concise_",tsv_files[n], sep=""), sep=""), sep="\t")
          
          ## remove objects to free up space
          rm(runner_raw_concise_data)
     }
          
          # # identify camera fov
          # ## plot trail outline
          # plot(clean_trail_data[,"X"], clean_trail_data[,"Y"], xlab="X", ylab="Y")
          # ## 4th camera
          # points(clean_trail_data[c(3,14,46,58),"X"],clean_trail_data[c(3,14,46,58),"Y"], pch=21, bg="blue")
          # polygon(clean_trail_data[c(14,3,46,58),"X"],clean_trail_data[c(14,3,46,58),"Y"],
          #         density=NA, col=scales::alpha("blue", alpha=0.3))
          # ## 3rd camera
          # points(clean_trail_data[c(25,49,71,79),"X"],clean_trail_data[c(25,49,71,79),"Y"], pch=21, bg="blue")
          # polygon(clean_trail_data[c(25,49,71,79),"X"],clean_trail_data[c(25,49,71,79),"Y"],
          #         density=NA, col=scales::alpha("blue", alpha=0.3))
          # ## 2nd camera
          # points(clean_trail_data[c(27,32,35,59),"X"],clean_trail_data[c(27,32,35,59),"Y"], pch=21, bg="blue")
          # polygon(clean_trail_data[c(32,27,35,59),"X"],clean_trail_data[c(32,27,35,59),"Y"],
          #         density=NA, col=scales::alpha("blue", alpha=0.3))
          # ## 1st camera
          # points(clean_trail_data[c(29,33,40,64),"X"],clean_trail_data[c(29,33,40,64),"Y"], pch=21, bg="blue")
          # polygon(clean_trail_data[c(33,40,29,64),"X"],clean_trail_data[c(33,40,29,64),"Y"],
          #         density=NA, col=scales::alpha("blue", alpha=0.3))
          # 
          # ## save camera fov in matrix
          # camera_fov = matrix(NA, nrow = 16, ncol=4); colnames(camera_fov) = c("X", "Y", "ID", "pos")
          # camera_fov[,"ID"] = sort(rep(1:4,4)); camera_fov[,"pos"] = rep(1:4,4)
          # camera_fov[,"X"] = c(clean_trail_data[c(64,29,33,40, 35,59,27,32, 49,71,25,79, 14,3,58,46),"X"])
          # camera_fov[,"Y"] = c(clean_trail_data[c(64,29,33,40, 35,59,27,32, 49,71,25,79, 14,3,58,46),"Y"])
          # 
          # ## export matrix
          # write.table(camera_fov, file=paste(dir_data,"camera_fov.txt", sep=""), sep="\t")
          # 
          # ## remove object
          # rm(camera_fov)
          
     ## extract first and last datapoint for each tracklet ####
          
          ## load data and make matrix
          concise_data = read.table(paste(dir_data,"concise_Quality_check0001_04102022.tsv", sep="")); concise_data = as.matrix(concise_data); rownames(concise_data) = NULL
          
          ## extract first and last item from each tracklet
          first_last_list = first.last.finder(concise_data)
          
          ## export short dataframes
          write.table(first_last_list[[1]], file=paste(dir_data,"Quality_check0001_04102022_stitching_data_first.txt", sep=""), sep="\t")
          write.table(first_last_list[[2]], file=paste(dir_data,"Quality_check0001_04102022_stitching_data_last.txt", sep=""), sep="\t")
          
     ## stitch tracklets together (5000 iterations < 3 mins) ####
          
          ## load data and make matrix for faster computation
          stitching_data_first = read.table(paste(dir_data,"Quality_check0001_04102022_stitching_data_first.txt", sep="")) 
          stitching_data_last = read.table(paste(dir_data,"Quality_check0001_04102022_stitching_data_last.txt", sep=""))
          stitching_data_first = as.matrix(stitching_data_first); rownames(stitching_data_first) = NULL
          stitching_data_last = as.matrix(stitching_data_last); rownames(stitching_data_last) = NULL
          
          ## list for track_IDs
          track_ID = vector(mode='list', length=max(unique(stitching_data_last[,"ID"])))
          
          ## sequential stitching
          seq_results = sequential.fragment.stitcher(u_first = stitching_data_first,
                                                     u_last = stitching_data_last,
                                                     u_track_ID = track_ID,
                                                     u_time_window_start = 1,
                                                     u_time_window_steps = 0.1,
                                                     u_x_diff_start = 1,
                                                     u_x_diff_steps = 0.1,
                                                     u_y_diff_start = 1,
                                                     u_y_diff_steps = 0.1,
                                                     iterations = 5000)
          
          ## extract tracks
          track_list = track.extractor(u_track_ID = seq_results[[3]])
          
          ## assemble full track data
          concise_data = read.table(paste(dir_data,"concise_Quality_check0001_04102022.tsv", sep=""));concise_data = as.matrix(concise_data); rownames(concise_data) = NULL
          track_data = track.data(u_track_ID = track_list, u_data = concise_data)
          
          ## remove tracks shorter than 10 minutes
          ## track_data = track_data[which(unlist(lapply(track_data,nrow)) > 6000)]
          
     ## testing ####
          
          ## Are there overlapping time frame values in the tracks?
          overlapping_time_test = function(u_track_data){
               
               if ((sum(sapply(u_track_data, function(x) length(unique(x[,"frame_number"])) == nrow(x))) != length(track_data)) == FALSE) {
                    print("TEST PASSED: There are no overlapping time frame values")
               } else {
                    print("TEST FAILED: There are overlapping time frame values")
                    
               }
          }
          overlapping_time_test(track_data)
          
          ## Is any tracklet used multiple times?
          tracklets_used_multiple_times_test = function(u_track_data) {
               if (sum(duplicated(unlist(lapply(u_track_data, function(x) unique(x[,"ID"]))))) == 0)
                    print("TEST PASSED: no tracklet was used in more than one track")
               else {
                    print("TEST FAILED: at least one tracklet was used on two tracks")
               }
               
          }
          tracklets_used_multiple_times_test(track_data)
          
          ## Any tracklets fell through the stitching process?
          unstitched_tracklets_test = function(u_track_data, u_raw_data){
               
               ## tracklets from tracks and raw data
               from_tracks = length(unlist(lapply(u_track_data, function(x) unique(x[,"ID"]))))
               from_fulldata = length(unique(u_raw_data[,"ID"]))
               
               if (length(from_tracks) == length(from_fulldata)){
                    print("TEST PASSED: all tracklets have been used in the stitching process")
               } else {
                    print("TEST FAILED: the number of tracklets in raw data does not match the number of tracklets post stitching")
               }
               
          }
          unstitched_tracklets_test(track_data, concise_data)
          
          ## Does the distance (normalized by time) between tracklets exceed 30?
          stitched_distance_test = function(u_track_data, max_stitch_distance) {
               
               distance_list = vector(mode='list', length=length(u_track_data))
               for (n in 1:length(u_track_data)) {
                    
                    ## run through IDs
                    runner_track = u_track_data[[n]]
                    
                    ## find first instance of tracklet
                    first_vec = vector()
                    for (m in unique(runner_track[,"ID"])) first_vec = c(first_vec, head(which(runner_track[,"ID"] == m),1))
                    
                    ## find last instance of tracklet
                    last_vec = vector()
                    for (k in unique(runner_track[,"ID"])) last_vec = c(last_vec, tail(which(runner_track[,"ID"] == k),1))
                    
                    ## calculate distance between last and first
                    distance_vec = vector()
                    for (o in 1:length(last_vec)) {
                         
                         ## end of track; start of track
                         runner_end = runner_track[last_vec[o],,drop=FALSE]; runner_start = runner_track[first_vec[o+1],,drop=FALSE]
                         
                         ## calculate distance: sqrt((x1-x2)^2+(y1-y2)^2)
                         raw_distance = round( sqrt((as.numeric(runner_start[,"X"])-as.numeric(runner_end[,"X"]))^2+(as.numeric(runner_start[,"Y"])-as.numeric(runner_end[,"Y"]))^2) ,3)
                         
                         ## calculate difference in time
                         time_dif = as.numeric(runner_start[,"frame_number"] - runner_end[,"frame_number"])
                         
                         ## normalize by time
                         distance_vec[o] = round(raw_distance/time_dif,3)
                    }
                    
                    ## collect distances in list
                    distance_list[[n]] = distance_vec
               }
          
               ## which tracks exceed the max_stitch_distance?
               flagged_tracks = lapply(distance_list, function(x) max(x, na.rm = TRUE) >= max_stitch_distance)
               flagged_tracks = which(flagged_tracks == TRUE)
               
               ## output
               if (length(flagged_tracks) == 0){
                    print(paste("TEST PASSED: no distance between stitched tracklets exceeds ",max_stitch_distance,"mm per time unit", sep=""))
               } else {
                    print(paste("TEST FAILED: these tracks", paste(flagged_tracks,collapse = ","), " exceed the distance of ",
                                max_stitch_distance, "mm per time unit between stitched tracklets. Check stitching paramters"))
               }
               
          }
          stitched_distance_test(track_data, 50)
          
          ## how long are the tracks
          plot(1:length(unlist(lapply(track_data, nrow))), sort(unlist(lapply(track_data, nrow))))
          
          # # plot tracks
          # for (n in 1:length(track_data)) {
          # 
          #      runner_track = track_data[[n]]
          #      palette(RColorBrewer::brewer.pal(12, "Set3"))
          # 
          #      ## create file for export
          #      png(filename=paste("Track_",n,".png",sep=""))
          # 
          #      plot(runner_track[,"X"], runner_track[,"Y"], col=runner_track[,"ID"],pch=16, cex=1)
          # 
          #      ## close graphic device
          #      dev.off()
          # }
          
## 2. DATA ANALYSIS (WAITING TIMES) ####
     ## how often are tracks near to each other? ####
          
          ## extract stitch events
          stitch_extractor = function(x) {
               
               stitching_events = which(diff(x[,"ID"]) != 0)
               first_track = x[stitching_events+1,, drop=FALSE]
               last_track = x[stitching_events,,drop=FALSE]     
               
               return_list = list(first_track, last_track)
               return(return_list)
          }
          stitch_events = lapply(track_data,function(x) stitch_extractor(x))
          
          ## make master matrix for last
          stitch_last = matrix(nrow=0, ncol=length(colnames(as.matrix(as.data.frame(stitch_events[[1]][2])))),dimnames = list(NULL,colnames(as.matrix(as.data.frame(stitch_events[[1]][2])))))
          for(n in 1:length(stitch_events)) stitch_last = rbind(stitch_last, as.data.frame(stitch_events[[n]][2]))
          
          ## make master matrix for first
          stitch_first = matrix(nrow=0, ncol=length(colnames(as.matrix(as.data.frame(stitch_events[[1]][2])))),dimnames = list(NULL,colnames(as.matrix(as.data.frame(stitch_events[[1]][2])))))
          for(n in 1:length(stitch_events)) stitch_first = rbind(stitch_first, as.matrix(as.data.frame(stitch_events[[n]][2])))
          
          ## workhorse
          summary_list = vector(mode="list", length=length(unique(stitch_last[,"ID_track"])))
          for (n in 1:length(unique(stitch_last[,"ID_track"]))) {
               
               ## progress bar
               print(round(n/length(unique(stitch_last[,"ID_track"])),4)*100)
               
               ## track IDs
               track_IDs = unique(stitch_last[,"ID_track"])
               
               ## extract first track
               runner_last = stitch_last[which(stitch_last[,"ID_track"] == track_IDs[n]),]
               
               ## find cases where x,y,t are similar between runner_last and all firsts
               similar_x = vector(mode="list", length(nrow(runner_last)))
               similar_y = vector(mode="list", length(nrow(runner_last)))
               similar_t = vector(mode="list", length(nrow(runner_last)))
               
               for(m in 1:nrow(runner_last)) similar_x[[m]] = which(lapply(stitch_first[,"X"],function(x) abs(as.numeric(x)-as.numeric(runner_last[m,"X"]))) <= 100)
               for(m in 1:nrow(runner_last)) similar_y[[m]] = which(lapply(stitch_first[,"Y"],function(x) abs(as.numeric(x)-as.numeric(runner_last[m,"Y"]))) <= 100)
               for(m in 1:nrow(runner_last)) similar_t[[m]] = which(lapply(stitch_first[,"frame_number"],function(x) abs(as.numeric(x)-as.numeric(runner_last[m,"frame_number"]))) <= 100)
               
               ## find intersection between x,y,t
               candidates = vector(mode="list", length(nrow(runner_last)))
               for(m in 1:nrow(runner_last)) candidates[[m]] = intersect(similar_x[[m]], intersect(similar_y[[m]], similar_t[[m]]))
               
               ## exclude correct stitching
               candidates_matrix = vector(mode="list", length(nrow(runner_last)))
               for(m in 1:nrow(runner_last)) candidates_matrix[[m]] = stitch_first[intersect(similar_x[[m]], intersect(similar_y[[m]], similar_t[[m]])),, drop=FALSE]
               summary_list[[n]] = which(lapply(lapply(candidates_matrix, function(x) x[,"ID_track"] != unique(runner_last[,"ID_track"])), sum) > 1)
               
          }
          
          tracks_without_conflict = which(unlist(lapply(summary_list,sum)) == 0)
          
          plot(track_data[[81]][,"X"], track_data[[81]][,"Y"])
          points(track_data[[10]][,"X"], track_data[[10]][,"Y"])
          
          
          
          
          stoppter
          
     ## extract details for video matching ####
          
          ## load camera_fov
          camera_fov = read.table(file = paste(dir_data, "camera_fov.txt", sep=""), sep = '\t')
          
          ## create fovs cutoffs
          ## format low x, high x, low y, high y
          camera_cutoff = vector(mode='list', length=4)
          camera_cutoff[[1]] = c(1300, 2750, 1800, 2250)
          camera_cutoff[[2]] = c(-1300, 100, -2350, -1900)
          camera_cutoff[[3]] = c(3500, 4500, 4350, 5400)
          camera_cutoff[[4]] = c(-1850, -400, 8020, 8410)
          
          ## checks if tracks are under camera, returns time spans in minutes
          under.camera = function(u_track_data, u_camera_cutoff) {
               
               ## output list
               final_output = matrix(nrow = length(u_track_data), ncol=9) 
               colnames(final_output) = c("track_ID", "c1_start", "c1_end", "c2_start", "c2_end",
                                                      "c3_start", "c3_end", "c4_start", "c4_end")
                                          
               for (n in 1:length(u_track_data)) {
                    
                    runner_track = u_track_data[[n]]
               
                    ## test for cameras
                    camera_one = which(runner_track[,"X"] >= u_camera_cutoff[[1]][1] & runner_track[,"X"] <= u_camera_cutoff[[1]][2] &
                                       runner_track[,"Y"] >= u_camera_cutoff[[1]][3] & runner_track[,"Y"] <= u_camera_cutoff[[1]][4])
                    camera_two = which(runner_track[,"X"] >= u_camera_cutoff[[2]][1] & runner_track[,"X"] <= u_camera_cutoff[[2]][2] &
                                       runner_track[,"Y"] >= u_camera_cutoff[[2]][3] & runner_track[,"Y"] <= u_camera_cutoff[[2]][4])
                    camera_three = which(runner_track[,"X"] >= u_camera_cutoff[[3]][1] & runner_track[,"X"] <= u_camera_cutoff[[3]][2] &
                                         runner_track[,"Y"] >= u_camera_cutoff[[3]][3] & runner_track[,"Y"] <= u_camera_cutoff[[3]][4])
                    camera_four = which(runner_track[,"X"] >= u_camera_cutoff[[4]][1] & runner_track[,"X"] <= u_camera_cutoff[[4]][2] &
                                        runner_track[,"Y"] >= u_camera_cutoff[[4]][3] & runner_track[,"Y"] <= u_camera_cutoff[[4]][4])
                    
                    ## extract times
                    camera_one_mins = round(range(runner_track[camera_one,"time_in_deciseconds"])/60,2)
                    camera_two_mins = round(range(runner_track[camera_two,"time_in_deciseconds"])/60,2)
                    camera_three_mins = round(range(runner_track[camera_three,"time_in_deciseconds"])/60,2)
                    camera_four_mins = round(range(runner_track[camera_four,"time_in_deciseconds"])/60,2)
               
                    ## fill output matrix
                    final_output[n,"track_ID"] = n
                    final_output[n,"c1_start"] = camera_one_mins[1]; final_output[n,"c1_end"] = camera_one_mins[2]
                    final_output[n,"c2_start"] = camera_two_mins[1]; final_output[n,"c2_end"] = camera_two_mins[2]
                    final_output[n,"c3_start"] = camera_three_mins[1]; final_output[n,"c3_end"] = camera_three_mins[2]
                    final_output[n,"c4_start"] = camera_four_mins[1]; final_output[n,"c4_end"] = camera_four_mins[2]
                    
               }
               
               ## replace Inf and -Inf by NA
               final_output[sapply(final_output, is.infinite)] <- NA
               
               ## return output
               return(final_output)
          }
          
          ## run function
          under_camera = under.camera(track_data, camera_cutoff)
          
          ## identify tracks that fall under camera
          number_of_cameras_per_track = lapply(under_camera, function(x) (sum(is.infinite(unlist(x)) == FALSE)/2))
          usable_tracks = which(number_of_cameras_per_track != 0)
          
          # print(paste(round(length(usable_tracks)/length(track_data),2)*100,"% of tracks fall under at least one camera"), sep="")
          
          ## create matrix to check stitching via cameras
          
     ## calculate waiting times #### 
          
          ## list to store waiting events
          waiting_data = vector(mode='list', length=length(track_data))
          
          ## extract waiting events
          for (n in 1:length(track_data)) {
               
               ## run through tracks
               runner_track = track_data[[n]]
               
               ## remove NAs from tracks by computing average of the two adjacent values
               runner_na = which(is.na(runner_track[,"distance"]) == TRUE)
               for (m in 1:(length(runner_na)-1)) runner_track[runner_na[m], "distance"] = mean(c(runner_track[runner_na[m]-1,"distance"], runner_track[runner_na[m]+1,"distance"]))
               
               ## smooth distances
               for (m in 1:4) runner_track[,"distance"] = runmed(runner_track[,"distance"],111)
               
               windows(50,15)
               plot(1:nrow(runner_track), log(runner_track[,"distance"]), type="l")
               plot(1:1000, log(runner_track[1:1000,"distance"]), type="l")
               abline(h=log(1.5))
               
               ## set cutoff at 1.5mm/ds which equals 15 mm/s which equals 30% of average walking speed (https://doi.org/10.1007/s000400300001)
               walking_track = which(runner_track[,"distance"] > 1.5)
               resting_track = which(runner_track[,"distance"] < 1.5)
               
               plot(runner_track[,"X"], runner_track[,"Y"], type="l")
               points(runner_track[resting_track,"X"], runner_track[resting_track,"Y"],col="blue")
               points(runner_track[walking_track,"X"], runner_track[walking_track,"Y"], type="l",col="red")
               
               ## extract separate resting events (allow movement for 100 seconds during events)
               sep_events = c(1,which(diff(runner_track[resting_track,"frame_number"]) > 100))
               
               ## extract length of each event
               sep_len = diff(sep_events)
               
               ## remove events shorter than 10 seconds
               waiting_data[[n]] = sep_len[which(sep_len > 100)]
          }
          
          
## 3. Video data ####
     ## load and handle data ####
     
     ## load raw data
     video_data = read.table(file = paste("C:/Users/timte/Desktop/Konstanz/Chain transport/Erik/master/Ant_Trail_20221007/master_file_3_20221007.txt"),sep = '\t', header=TRUE, fill=FALSE)

     ## transform master file into convenient table for matching to tracks
          
          ## create list for split dataframes
          runner_dfs = vector(mode="list", length=3)
          
          ## extract brakes
          breaks = which(video_data[,10] == "BREAK")
          
          ## split dataframes
          if (length(breaks) == 0) { runner_dfs[[1]] = video_data
          } else if (length(breaks) == 1){ runner_dfs[[1]] = video_data[1:breaks[1],]; runner_dfs[[2]] = video_data[breaks[1]:nrow(video_data),]
          } else if (length(breaks) == 2){ runner_dfs[[1]] = video_data[1:breaks[1],]; runner_dfs[[2]] = video_data[breaks[1]:breaks[2],]; runner_dfs[[3]] = video_data[breaks[2]:nrow(video_data),]}
          
          ## reset time col in the split dfs
          for (m in 1:(length(breaks)+1)) runner_dfs[[m]][,"Time"] = 1:nrow(runner_dfs[[m]])
          
          ## find shape in each col
          shape_by_col = vector(mode="list", length=length(breaks)+1)
          for (m in 1:(length(breaks)+1)) shape_by_col[[m]] = apply(runner_dfs[[m]][,2:ncol(runner_dfs[[m]])], 2, function(x) which(x %in% c("STAR", "HEART", "CIRCLE") == TRUE))
          for (m in 1:(length(breaks)+1)) shape_by_col[[m]] = as.numeric(which(shape_by_col[[m]] > 0))+1
          
          ## create output 
          output_list_df = vector(mode="list", length=length(breaks)+1)
          
          ## extract shape and store time and shape in df
          for (k in 1:(length(breaks)+1)) {
               
               ## extract cols with a shape
               runner = runner_dfs[[k]][,c(1,shape_by_col[[k]])]
               
               ## extract entry time
               runner_time =  runner[as.numeric(apply(runner[,2:ncol(runner)], 2, function(x) which(x == "ENTER_RIGHT"))),"Time"]
               
               ## extract shape
               runner_shape = runner[apply(runner[,2:ncol(runner)], 2, function(x) which(x %in% c("STAR", "HEART", "CIRCLE"))),2:ncol(runner)]
               runner_shape = as.vector(unlist(runner_shape)[which(unlist(runner_shape) != "")])
               runner_shape = runner_shape[which(runner_shape %in% c("STAR", "HEART", "CIRCLE"))]
               
               ## combine in df
               output_list_df[[k]] = data.frame("time" = runner_time,
                                                "shape" = runner_shape)
          }
     

     
     
## 4. ARCHIVE ####
     ## plotter ####
     
     # plot tracks
     for (n in 1:length(track_data)) {
          
          runner_track = track_data[[n]]
          
          ## create file for export
          png(filename=paste("Distances in Track_",n,".png",sep=""), width=4000, height=1000)
          
          plot(1:nrow(runner_track), log(runner_track[,"distance"]), type="l")
          
          for (k in which(is.na(runner_track[,"distance"]) == TRUE)) abline(v=k,col="red")
          
          ## close graphic device
          dev.off()
     }
     
     ## stitching example that works ####
     
     ## version that works for the example
     fragment.stitcher = function (u_first, u_last, u_time_window, u_x_diff, u_y_diff) {
          
          # for testing
          # u_first = fake_data_first; u_last = fake_data_last; u_time_window = 5; u_x_diff = 3; u_y_diff = 3; n=1
          
          ## list for track_IDs
          track_ID = vector(mode='list', length=length(unique(u_last[,"ID"])))
          
          ## for loop (runs through all endpoints IDs)
          ## for (n in 1:length(unique(u_last[,"ID"]))){
          
          ## variables for the while loop
          proceed = TRUE
          n=1
          
          while (proceed == TRUE) {
               
               ## take endpoint
               runner_endpoint = u_last[which(u_last[,"ID"] == n),,drop=FALSE]
               
               ## select all startpoints, that the come after runner_endpoint and are within u_time_window
               runner_candidates_time = u_first[which(u_first[,"frame_number"] > as.numeric(runner_endpoint[,"frame_number"]) & 
                                                           u_first[,"frame_number"] <= (as.numeric(runner_endpoint[,"frame_number"])+u_time_window)),, drop=FALSE]
               
               ## select all rows differ less than u_x_diff units on X
               runner_candidate_time_x = runner_candidates_time[which((abs(runner_candidates_time[,"X"] - as.numeric(runner_endpoint[,"X"]))) <= u_x_diff),, drop=FALSE]
               
               ## select all rows differ less than u_y_diff units on Y
               runner_candidate_time_x_y = runner_candidate_time_x[which((abs(runner_candidate_time_x[,"Y"] - as.numeric(runner_endpoint[,"Y"]))) <= u_y_diff),, drop=FALSE]
               
               ## if one tracklet left: stitch
               if (nrow(runner_candidate_time_x_y) == 1) {
                    
                    ## update track_ID
                    track_ID[[runner_endpoint[,"ID"]]] = paste(track_ID[[n]],runner_endpoint[,"ID"],runner_candidate_time_x_y[,"ID"],sep=",")
                    
                    ## remove endpoint of runner_endpoint
                    u_last = u_last[-which(u_last[,"ID"] == runner_endpoint[,"ID"]),, drop=FALSE]
                    
                    ## update new endpoint ID
                    u_last[which(u_last[,"ID"] == runner_candidate_time_x_y[,"ID"]),"ID"] = runner_endpoint[,"ID"]
                    
                    ## remove now sitched startpoint runner_candidate_time_x_y
                    u_first = u_first[-which(u_first[,"ID"] == runner_candidate_time_x_y[,"ID"]),, drop=FALSE]
                    
               } else {
                    
                    ## go to next tracklet
                    n = n + 1
                    
                    ## stop if there is no more next tracklet
                    proceed = ifelse(n > length(u_last[,"ID"]), FALSE, TRUE)
               }
               
          }
          
          ## make return list
          return_list = list(u_first, u_last, track_ID)     
          
          ## function output
          return(return_list)
     }
     
     ## create scenario
     fake_data_first = stitching_data_first
     fake_data_last = stitching_data_last
     
     fake_data_first[1,] = c(1,1,2,2,100,1)
     fake_data_first[2,] = c(3,3,10,8,100,2)
     fake_data_first[3,] = c(5,5,19,13,100,3)
     fake_data_first[4,] = c(7,7,24,20,100,4)
     fake_data_first[5,] = c(9,9,23,11,100,5)
     fake_data_first[6,] = c(11,11,16,15,100,6)
     fake_data_first = fake_data_first[-c(7:nrow(fake_data_first)),]
     
     fake_data_last[1,] = c(2,2,9,6,100,1)
     fake_data_last[2,] = c(4,4,18,11,100,2)
     fake_data_last[3,] = c(6,6,23,18,100,3)
     fake_data_last[4,] = c(8,8,24,12,100,4)
     fake_data_last[5,] = c(10,10,16,15,100,5)
     fake_data_last[6,] = c(12,12,11,10,100,6)
     fake_data_last = fake_data_last[-c(7:nrow(fake_data_last)),]
     
     plot(fake_data_first[,"X"], fake_data_first[,"Y"], col="dodgerblue", pch=16)
     points(fake_data_last[,"X"], fake_data_last[,"Y"], col="darkgoldenrod", pch=16)
     for (n in 1:6) lines(c(fake_data_first[,"X"][n], fake_data_last[,"X"][n]),c(fake_data_first[,"Y"][n], fake_data_last[,"Y"][n]), lty=2)
     legend("topleft", inset=0.02, legend=c("startpoint", "endpoint"), pch=16,
            col = c("dodgerblue", "darkgoldenrod"))
     
     ## list for track_IDs
     track_ID = vector(mode='list', length=length(unique(fake_data_last[,"ID"])))
     
     kek = sequential.fragment.stitcher(u_first = fake_data_first,
                                        u_last = fake_data_last,
                                        u_track_ID = track_ID,
                                        u_time_window_start = 1,
                                        u_time_window_steps = 0.3,
                                        u_x_diff_start = 1,
                                        u_x_diff_steps = 0.4,
                                        u_y_diff_start = 1,
                                        u_y_diff_steps = 0.4,
                                        iterations = 10)
                         
     
     ## old test data ####
     
     ## make a fake list
     fake_list = vector(mode='list', length=5)
     fake_list[[1]] = NULL
     fake_list[[2]] = c(2,4)
     fake_list[[3]] = NULL
     fake_list[[4]] = c(4,5,6)
     fake_list[[5]] = NULL
     fake_list[[6]] = c(6,7)
     fake_list[[7]] = NULL
     
     ## generate new test data
     fake_data_first = stitching_data_first
     fake_data_last = stitching_data_last
     
     fake_data_first[1,] = c(1,1,2,2,100,1)
     fake_data_first[2,] = c(3,3,10,8,100,2)
     fake_data_first[3,] = c(5,5,19,13,100,3)
     fake_data_first[4,] = c(7,7,24,20,100,4)
     fake_data_first[5,] = c(9,9,23,11,100,5)
     fake_data_first[6,] = c(11,11,16,15,100,6)
     fake_data_first = fake_data_first[-c(7:nrow(fake_data_first)),]
     
     fake_data_last[1,] = c(2,2,9,6,100,1)
     fake_data_last[2,] = c(4,4,18,11,100,2)
     fake_data_last[3,] = c(6,6,23,18,100,3)
     fake_data_last[4,] = c(8,8,24,12,100,4)
     fake_data_last[5,] = c(10,10,16,15,100,5)
     fake_data_last[6,] = c(12,12,11,10,100,6)
     fake_data_last = fake_data_last[-c(7:nrow(fake_data_last)),]
     
     plot(fake_data_first[,"X"], fake_data_first[,"Y"], col="dodgerblue", pch=16,xlim=c(0,25), ylim=c(0,20))
     points(fake_data_last[,"X"], fake_data_last[,"Y"], col="darkgoldenrod", pch=16)
     for (n in 1:6) lines(c(fake_data_first[,"X"][n], fake_data_last[,"X"][n]),c(fake_data_first[,"Y"][n], fake_data_last[,"Y"][n]), lty=2)
     legend("topleft", legend=c("startpoint", "endpoint"), pch=16,
            col = c("dodgerblue", "darkgoldenrod"))
     for (n in 1:25) abline(v=n, col=scales::alpha("black", 0.1)); for (n in 1:20) abline(h=n, col=scales::alpha("black", 0.1))
     
     
     track_ID = vector(mode='list', length=length(unique(fake_data_last[,"ID"])))
     kek = fragment.stitcher(u_first = fake_data_first,
                             u_last = fake_data_last,
                             u_track_ID = track_ID,
                             u_time_window = 1,
                             u_x_diff = 2,
                             u_y_diff = 2)
     
     
     
     
     
     fake_results = sequential.fragment.stitcher(u_first = fake_data_first,
                                                 u_last = fake_data_last,
                                                 u_track_ID = track_ID,
                                                 u_time_window_start = 1,
                                                 u_time_window_steps = 0,
                                                 u_x_diff_start = 0.5,
                                                 u_x_diff_steps = 0.5,
                                                 u_y_diff_start = 0.5,
                                                 u_y_diff_steps = 0.5,
                                                 iterations = 3)
     
     iter_one = fragment.stitcher(u_first = fake_data_first,
                                  u_last = fake_data_last,
                                  u_track_ID = track_ID,
                                  u_time_window = 1,
                                  u_x_diff = 1,
                                  u_y_diff = 1)
     
     iter_two = fragment.stitcher(u_first = iter_one[[1]],
                                  u_last = iter_one[[2]],
                                  u_track_ID = iter_one[[3]],
                                  u_time_window = 1,
                                  u_x_diff = 1.5,
                                  u_y_diff = 1.5)
     
     iter_three = fragment.stitcher(u_first = iter_two[[1]],
                                    u_last = iter_two[[2]],
                                    u_track_ID = iter_two[[3]],
                                    u_time_window = 1,
                                    u_x_diff = 2,
                                    u_y_diff = 2)
     
     
     
     ## random old stuff ####
     
     still_works = sequential.fragment.stitcher(u_first = stitching_data_first,
                                                u_last = stitching_data_last,
                                                u_time_window_start = 1,
                                                u_time_window_steps = 1,
                                                u_x_diff_start = 1,
                                                u_x_diff_steps = 0.2,
                                                u_y_diff_start = 1,
                                                u_y_diff_steps = 0.2,
                                                iterations = 10)
     
     ## create test data
     tracklet_1 = data.frame("frame_number"=1:10,   "time_in_deciseconds"=(1:10)/10,   "X"=101:110,   "Y"=101:110,   "Z"=101:110,   "ID"= rep(1,10)  ); tracklet_1 = as.matrix(tracklet_1)
     tracklet_2 = data.frame("frame_number"=12:20,  "time_in_deciseconds"=(12:20)/10,  "X"=112:120,   "Y"=112:120,   "Z"=112:120,   "ID"= rep(2,9)   ); tracklet_2 = as.matrix(tracklet_2)
     tracklet_3 = data.frame("frame_number"=25:30,  "time_in_deciseconds"=(25:30)/10,  "X"=115:120,   "Y"=115:120,   "Z"=115:120,   "ID"= rep(3,6)   ); tracklet_3 = as.matrix(tracklet_3)
     
     ## wrong place, right time
     tracklet_4 = data.frame("frame_number"=32:40,  "time_in_deciseconds"=(32:40)/10,  "X"=32:40,     "Y"=32:40,     "Z"=32:40,     "ID"= rep(4,9)   ); tracklet_4 = as.matrix(tracklet_4)
     
     ## right place, wrong time
     tracklet_5 = data.frame("frame_number"=1:11,  "time_in_deciseconds"=(1:11)/10,  "X"=120:130,     "Y"=120:130,     "Z"=120:130,     "ID"= rep(5,11)   ); tracklet_5 = as.matrix(tracklet_5)
     
     ## combine to test data
     test_pretty_data = rbind(tracklet_1, tracklet_2, tracklet_3, tracklet_4, tracklet_5)
     
     ## generate first_and_last
     test_first_last = first.last.finder(test_pretty_data)
     
     
     
     ## test stitching
     kek = sequential.fragment.stitcher(u_first = test_first_last[[1]],
                                        u_last = test_first_last[[2]],
                                        u_time_window_start = 1,
                                        u_time_window_steps = 1,
                                        u_x_diff_start = 1,
                                        u_x_diff_steps = 1,
                                        u_y_diff_start = 1,
                                        u_y_diff_steps = 1,
                                        iterations = 10)
     
     
     
     
     iter_1 = fragment.stitcher(stitching_data_first, stitching_data_last, 1,1,1)
     
     
     
     iter_2 = fragment.stitcher(iter_1[[1]], iter_1[[2]], 2,2,2)
     iter_3 = fragment.stitcher(iter_2[[1]], iter_2[[2]], 3,3,3)
     iter_4 = fragment.stitcher(iter_3[[1]], iter_3[[2]], 4,4,4)
     iter_5 = fragment.stitcher(iter_4[[1]], iter_4[[2]], 5,5,5)
     iter_6 = fragment.stitcher(iter_5[[1]], iter_5[[2]], 6,6,6)
     iter_7 = fragment.stitcher(iter_6[[1]], iter_6[[2]], 7,7,7)
     iter_8 = fragment.stitcher(iter_7[[1]], iter_7[[2]], 8,8,8)
     iter_9 = fragment.stitcher(iter_8[[1]], iter_8[[2]], 9,9,9)
     iter_10 = fragment.stitcher(iter_9[[1]], iter_9[[2]], 10,10,10)
     iter_11 = fragment.stitcher(iter_10[[1]], iter_10[[2]], 11,11,11)
     iter_12 = fragment.stitcher(iter_11[[1]], iter_11[[2]], 12,12,12)
     iter_13 = fragment.stitcher(iter_12[[1]], iter_12[[2]], 13,13,13)
     iter_14 = fragment.stitcher(iter_13[[1]], iter_13[[2]], 14,14,14)
     iter_15 = fragment.stitcher(iter_14[[1]], iter_14[[2]], 15,15,15)
     iter_16 = fragment.stitcher(iter_15[[1]], iter_15[[2]], 16,16,16)
     
     
     iter_16
     
     
     
     ## stitching of tracklets within a track must be chronological
     
     sub_track_first_tracklet = lapply(track_data, function(x) first.last.finder(x)[[1]])
     sub_track_last_tracklet = lapply(track_data, function(x) first.last.finder(x)[[2]])
     
     first_is_increasing = lapply(sub_track_first_tracklet, function(x) !is.unsorted(x[,"time_in_deciseconds"]))
     last_is_increasing = lapply(sub_track_last_tracklet, function(x) !is.unsorted(x[,"time_in_deciseconds"]))
     
     which(first_is_increasing==FALSE)
     which(last_is_increasing==FALSE)
     
     first_tn = sub_track_first_tracklet[[29]]
     runner_ID_subset= sub_track_last_tracklet[[29]]
     
     head(first_tn,5)
     head(last_tn,5)
     
     vec = vector()
     for (n in 1:nrow(first_tn)-1) vec[n] = last_tn[n,"time_in_deciseconds"] < first_tn[n+1,"time_in_deciseconds"]
     
     
     ## changes the format of the raw data into something more concise (all in dataframes) ####
     concise.dataframe = function(input_df) {
          
          ## labels the cols
          m=0
          names(input_df)[1] = "frame_number"
          names(input_df)[2] = "time_in_deciseconds"
          for(n in seq(from=3, to=(ncol(input_df)-2), by=3)) {
               
               ## counter to label the markers
               m = m + 1
               
               ## name each X,y & Z col with the counter m
               names(input_df)[n:(n+2)] = c(paste("X",m,sep=""), paste("Y",m,sep=""), paste("Z",m,sep=""))
          }
          
          ## create the new dataframe skeleton
          concise_df = data.frame("frame_number" = double(),
                                  "time_in_deciseconds" = double(),
                                  "X" = double(),
                                  "Y" = double(),
                                  "Z" = double(),
                                  "ID" = double())
          
          ## create ID counter
          ID_counter = 0
          
          ## reorganize dataframe
          for (n in seq(from=3, to=(ncol(input_df)-2), by=3)){
               
               ## add to ID counter
               ID_counter = ID_counter + 1
               
               ## extract coordinates for one marker
               runner_df = input_df[,c(1:2,n:(n+2))]
               
               ## find rows without NAs 
               good_rows = which(!is.na(runner_df[,3]))
               
               ## extract rows without NAs
               good_subset = runner_df[good_rows,]
               
               ## generate dataframe
               runner_merger = data.frame("frame_number" = good_subset$frame_number,
                                          "time_in_deciseconds" = good_subset$time_in_deciseconds,
                                          "X" = good_subset[,3],
                                          "Y" = good_subset[,4],
                                          "Z" = good_subset[,5],
                                          "ID" = rep(ID_counter,nrow(good_subset)))
               
               ## bind dataframe
               concise_df = rbind(concise_df, runner_merger)
          }
          
          return(concise_df)
          
     } 
     ## plot candidate tracklets ~ 1 minute ####
          
          
          ## extract full data for each tracklet and create plot
          for (n in 1:length(original_tracklets_for_each_stitched_tracklets)) {
               
               ## get ID
               runner_ID = names(original_tracklets_for_each_stitched_tracklets)[n]
               original_IDs = original_tracklets_for_each_stitched_tracklets[[n]]
               
               ## create plot with tracklet
               plotter_data = pretty_data[pretty_data[,"ID"] %in% original_IDs, ]
               
               ## sort by frame_number
               ## plotter_data[order(plotter_data[,"frame_number"],decreasing=FALSE),]
               
               ## create filename for export
               runner_filename = paste(runner_ID,".png", sep="")
               
               ## create file for export
               png(filename=runner_filename)
               
               ## create plot
               plot(plotter_data[,"X"], plotter_data[,"Y"],
                    col=plotter_data[,"ID"],
                    xlim=c(-3800,4600), ylim=c(-2200, 8300),
                    main=paste("Made from ",length(original_IDs)," tracklets", sep=""))
               
               ## close graphic device
               dev.off()
               
          }
          
     ## extract pretty tracks for TerraX ####
          
          ## extract data for three full tracklets (original IDs: 3569, 3260, 2407)
          terrax_one = pretty_data[pretty_data[,"ID"] %in% original_tracklets_for_each_stitched_tracklets$"3569", ]
          terrax_two = pretty_data[pretty_data[,"ID"] %in% original_tracklets_for_each_stitched_tracklets$"3260", ]
          terrax_three = pretty_data[pretty_data[,"ID"] %in% original_tracklets_for_each_stitched_tracklets$"2407", ]
          
          ## create a function
          distance.calculator = function(u_data) {
               
               ## full distance vector
               distance_vector = rep(NA, nrow(u_data))
               
               for (n in 1:length(unique(u_data[,"ID"]))) {
                    
                    ## subset the data by ID
                    runner_subset = u_data[which(u_data[,"ID"] == unique(u_data[,"ID"])[n]),]
                    
                    ## calculate distance
                    
                    ## vector to store distance
                    runner_vector = rep(NA, nrow(runner_subset))
                    
                    ## makes the calculations
                    for (m in 1:(nrow(runner_subset)-1)){
                         
                         ## get coordinates
                         x1 = as.numeric(runner_subset[m,"X"]);   y1 = as.numeric(runner_subset[m,"Y"])
                         x2 = as.numeric(runner_subset[m+1,"X"]); y2 = as.numeric(runner_subset[m+1,"Y"])
                         
                         ## calculate distance
                         runner_vector[m] = round( sqrt((x1-x2)^2+(y1-y2)^2) ,3)
                         
                    }
                    
                    ## add to final vector
                    distance_vector[which(u_data[,"ID"] == unique(u_data[,"ID"])[n])] = runner_vector
               }
               
               return(distance_vector)
          }
          
          ## get distance for all the tracks
          terrax_one = cbind(terrax_one, distance.calculator(terrax_one))
          terrax_two = cbind(terrax_two, distance.calculator(terrax_two))
          terrax_three = cbind(terrax_three, distance.calculator(terrax_three))
          
          ## remake ID
          terrax_one[,"ID"] = 1; terrax_two[,"ID"] = 2; terrax_three[,"ID"] = 3
          
          ## rename distance col
          colnames(terrax_one)[7] = "distance"; colnames(terrax_two)[7] = "distance"; colnames(terrax_three)[7] = "distance"
          
          ## combine data
          terrax_data = rbind(terrax_one, terrax_two, terrax_three)
          
          ## make look nice
          terrax_data = terrax_data[,colnames(terrax_data)!="frame_number"]; terrax_data = terrax_data[,colnames(terrax_data)!="Z"]
          colnames(terrax_data)[1] = "time_in_deciseconds_post_experiment_start"
          
          ## write.table(terrax_data,file=paste(raw_data_location,"TerraX_data.txt"), sep="\t")
          
          ## make plot     
          
          ## make palette
          cont_cols = brewer.pal(9, "YlOrRd")
          
          ## which track shall be drawn? (options are 1, 2 and 3)
          track_to_draw = 1
          
          ## draw plot
          plot(terrax_data[which(terrax_data[,"ID"] == track_to_draw),][,"X"],terrax_data[which(terrax_data[,"ID"] == track_to_draw),][,"Y"],
               col=ifelse(terrax_data[,"distance"] == 0, cont_cols[3],
                          ifelse(terrax_data[,"distance"] <= 1, cont_cols[4],
                                 ifelse(terrax_data[,"distance"] <= 5, cont_cols[7],
                                        ifelse(terrax_data[,"distance"] <= 10, cont_cols[9], cont_cols[9]))))
          )
          
          ## add legend
          legend("bottomright", inset=0.001, 
                 pch=16,
                 col = c(cont_cols[3],cont_cols[4],cont_cols[7],cont_cols[9],cont_cols[9]),
                 legend = c("0","<=1","<=5","<=10",">10"),
                 title="Speed: mm/deciseconds")
          
          
          
          ## stopper
     ## show trail outline ####
          
          ## load trail outline (need only first row)
          raw_trail_outline_data = read.table(file = paste(dir_data, "Trail_outline.tsv", sep=""), sep = '\t', header=FALSE, skip=12, fill=FALSE, nrows=2)
          
          ## transform to concise format
          clean_trail_data = m.concise.dataframe(raw_trail_outline_data)
          clean_trail_data = clean_trail_data[which(clean_trail_data[,"frame_number"] == 1),]
          
          
          ## plot trail outline
          plot(clean_trail_data[,"X"], clean_trail_data[,"Y"])
     
          ## 4th camera     
          points(clean_trail_data[c(3,14,46,58),"X"],clean_trail_data[c(3,14,46,58),"Y"], pch=21, bg="blue")
          polygon(clean_trail_data[c(14,3,46,58),"X"],clean_trail_data[c(14,3,46,58),"Y"],
                  density=NA, col=scales::alpha("blue", alpha=0.3))
          
          ## 3rd camera
          points(clean_trail_data[c(25,49,71,80),"X"],clean_trail_data[c(25,49,71,80),"Y"], pch=21, bg="blue")
          polygon(clean_trail_data[c(25,49,71,80),"X"],clean_trail_data[c(25,49,71,80),"Y"],
                  density=NA, col=scales::alpha("blue", alpha=0.3))
          
          ## 2nd camera
          points(clean_trail_data[c(27,32,35,59),"X"],clean_trail_data[c(27,32,35,59),"Y"], pch=21, bg="blue")
          polygon(clean_trail_data[c(32,27,35,59),"X"],clean_trail_data[c(32,27,35,59),"Y"],
                  density=NA, col=scales::alpha("blue", alpha=0.3))
          
          ## 1st camera
          points(clean_trail_data[c(29,33,40,64),"X"],clean_trail_data[c(29,33,40,64),"Y"], pch=21, bg="blue")
          polygon(clean_trail_data[c(33,40,29,64),"X"],clean_trail_data[c(33,40,29,64),"Y"],
                  density=NA, col=scales::alpha("blue", alpha=0.3))
          
          ## add text
          text(-3300,6500, "Colony")
          text(3300,2300, "Food")
          which(clean_trail_data[,"Z"] > 2000)
          
          
     ## fragment.stitcher tester ####
          
          
          # ## dummy data
          # tester_first = data.frame("frame_number" = c(1,7,30),
          #                           "time_in_deciseconds" = rep(NA,3),
          #                           "X" = c(10,16,201),
          #                           "Y" = c(1,6,312), 
          #                           "Z" = rep(NA,3),
          #                           "ID" = 1:3);tester_first = as.matrix(tester_first)
          # 
          # tester_last = data.frame("frame_number" = c(5,30,43),
          #                          "time_in_deciseconds" = rep(NA,3),
          #                          "X" = c(15,20,213),
          #                          "Y" = c(6,4,304), 
          #                          "Z" = rep(NA,3),
          #                          "ID" = 1:3);tester_last = as.matrix(tester_last)
          # 
          # ## to test my function manually 
          # u_first = tester_first; u_last = tester_last; u_time_window = 10; u_x_diff = 5
          # 
          ## to test my function
          # function_test = fragment.stitcher(u_first = tester_first, u_last = tester_last,
          #                                  u_time_window = 10, u_x_diff = 5)
          
          
     ## backup - before I tried to make it keep all tracklets under given ID ####
          
          fragment.stitcher = function (u_first, u_last, u_time_window, u_x_diff) {
               
               ## STEP 1. FIND MATCHES
               
               ## output vector
               matches_list = vector(mode='list', length=nrow(u_first))
               
               ## for loop
               for (n in 1:length(u_last[,"ID"])){
                    
                    ## take end of tracklet one by one
                    runner_row_last = u_last[n,, drop=FALSE]
                    
                    ## select all rows in u_first that are later than runner_row_last
                    runner_data = u_first[which(u_first[,"frame_number"] >= as.numeric(runner_row_last[,"frame_number"])),,drop=FALSE]
                    
                    ## select all rows in runner_data that come 100 frames after runner_row_last
                    runner_candidate_rows = runner_data[which(abs(runner_data[,"frame_number"] - as.numeric(runner_row_last[,"frame_number"])) <= u_time_window), , drop=FALSE]
                    
                    ## select all rows differ less than u_x_diff units on X
                    runner_candidate_rows = runner_candidate_rows[which((abs(runner_candidate_rows[,"X"] - as.numeric(runner_row_last[,"X"]))) <= u_x_diff),  , drop=FALSE]
                    
                    ## export candidate IDs
                    matches_list[[as.numeric(runner_row_last[,"ID"])]] = as.numeric(runner_candidate_rows[,"ID"])
               }     
               
               ## STEP 2. UPDATE DATA MATRICES WITH STITCHES
               
               ## extract the number of matches for each tracklet
               tracklet_ID_with_one_match = which(lapply(matches_list, function (x) length(x))==1)
               
               ## run through matches and update data
               for (n in tracklet_ID_with_one_match) {
                    
                    ## update label of match in first to n
                    u_first[which(u_first[,"ID"]==matches_list[[n]]), "ID"] = n
                    
                    ## update label of match in last to n
                    u_last[which(u_last[,"ID"]==matches_list[[n]]), "ID"] = n
                    
               }
               
               ## make return list
               return_list = list(u_first, u_last)
               
               ## function output
               return(return_list)
               
          }
          ## produces a first- and last-data matrix containing stitched tracklets
     ## old stitching #### 
          
          
          ## does the actual stitching and returns a stitched-dataframe 
          fragment.stitcher = function (u_data, matches_list) {
               
               ## use function output to combine fragments
               candidate_number = unlist(lapply(matches_list, function(x) length(x)))
               
               ## fragments which have only one candidate
               one_candidate_fragments = which(candidate_number == 1)
               
               ## make vector for remaining candidates
               unique(u_data[,"ID"])
               
               ## combine fragments with one candidate
               for (n in 1:length(one_candidate_fragments)) {
                    
                    ## take first fragment
                    runner_fragment = one_candidate_fragments[n]
                    
                    ## identify fragment that follows runner_fragment
                    runner_extender = matches_list[[runner_fragment]]
                    
                    ## change name of the extended fragment
                    u_data[which(u_data[,"ID"] == runner_extender),"ID"] = runner_fragment
               }
               return(u_data)
          }
          ## produces a "pretty_data" object with updated tracklet-IDs after stitching
          
          ## updates the first-last-matrices with the stitches 
          first.last.extractor = function (u_first, u_last) {
               
               
               
               
               
               
               
               ## store rownumbers
               first_index = rep(NA, length(unique(u_data[,"ID"])))
               last_index = rep(NA, length(unique(u_data[,"ID"])))     
               
               ## identify first and last datapoint of each tracklet
               for (n in 1:length(unique(u_data[,"ID"]))) {
                    
                    ## run through IDs
                    runner_ID = unique(u_data[,"ID"])[n]
                    
                    ## take last and first marker for each ID
                    first_index[n] = head(which(u_data[,"ID"] == runner_ID),1)
                    last_index[n] = tail(which(u_data[,"ID"] == runner_ID),1)
               }
               
               ## extract first and last datapoints
               stitching_data_first = pretty_data[first_index,]
               stitching_data_last = pretty_data[last_index,]
               
               ## sort matrix by frame number
               stitching_data_first = stitching_data_first[order(stitching_data_first[,1],decreasing = FALSE),]
               stitching_data_last = stitching_data_last[order(stitching_data_last[,1],decreasing = FALSE),]
               
               ## create list for return
               return_list = list(stitching_data_first, stitching_data_last)
               
               ## function ouput
               return(return_list)
          }
          ## produces a list with two matrices, the first containing the first and second the last datapoint for each tracklet
     ## old marker finding ####
          
          # ## remove tracks with clearly wrong Z values (i.e. 300 above or below Z=0)
          # ## removes 89 tracklets
          # 
          #    ## prepare output matrix
          #    extreme_z_values = matrix(nrow=0, ncol=3)
          #    colnames(extreme_z_values) = c("min","max","ID")
          #    
          #    ## find min and max Z values
          #    for (n in 1:length(unique(pretty_data[,"ID"]))){
          #         
          #         ## run through the ID vector
          #         runner_ID = unique(pretty_data[,"ID"])[n]
          #         
          #         ## extract values for ID
          #         runner_subset = pretty_data[which(pretty_data[,"ID"] == n),]
          #         
          #         ## get range for Z values
          #         runner_range = range(runner_subset[,"Z"])
          #         
          #         ## create runner matrix
          #         runner_matrix = matrix(c(runner_range[1], runner_range[2], runner_ID), ncol=3)
          #         
          #         ## bind matrix
          #         extreme_z_values = rbind(extreme_z_values, runner_matrix)
          #    }
          #    
          #    ## identify tracklets with with Z values above or below 1000mm
          #    to_remove_below_z = extreme_z_values[,"ID"][which(extreme_z_values[,"min"] <=-300)]
          #    to_remove_above_z = extreme_z_values[,"ID"][which(extreme_z_values[,"max"] >= 300)]
          #    to_remove_z = unique(c(to_remove_below_z, to_remove_above_z))
          #    
          #    ## remove extreme Z value tracklets from pretty data
          #    for (n in to_remove_z) pretty_data = pretty_data[-c(which(pretty_data[,"ID"] == n)),]
          #    
          #    ## remove Z col
          #    pretty_data = pretty_data[,colnames(pretty_data)!="Z"]
          #    
          # ## remove very short tracklets (i.e. less than 5 seconds)
          # ## removes 2783 tracklets
          #    
          #    ## get length of all tracklets
          #    tracklet_length_vector = vector()
          #    for (n in 1:length(unique(pretty_data[,"ID"]))) tracklet_length_vector[n] = length(which(pretty_data[,"ID"] == unique(pretty_data[,"ID"])[n]))
          #    
          #    ## find matching tracklet ID
          #    tracklets_to_remove_less_than_5_seconds = unique(pretty_data[,"ID"])[which(tracklet_length_vector <= 49)]
          #    tracklets_to_keep_more_than_5_seconds = setdiff(pretty_data[,"ID"], tracklets_to_remove_less_than_5_seconds)
          #    
          #    ## remove tracklets 
          #    ## for (n in tracklets_to_remove_less_than_5_seconds) pretty_data = pretty_data[-c(which(pretty_data[,"ID"] == n)),]
          #    pretty_data = pretty_data[pretty_data[,"ID"] %in% tracklets_to_keep_more_than_5_seconds,]
          
          ## remove tracklets that don't move more than 10cm in total (i.e. abs(max(x/y)-min(x/y)) <= 100)
          ## removes 1959
          
          # ## difference x values
          # x_y_diff = data.frame("x_diff" = double(),
          #                       "y_diff" = double(),
          #                       "ID" = double())
          # 
          # ## identify the change between max and min X value
          # for(n in 1:length(unique(pretty_data$ID))) {
          #      
          #      ## run through all IDs
          #      runner_ID = unique(pretty_data$ID)[n]
          #      
          #      ## calculate diff for X and Y values and store in dataframe
          #      runner_df = data.frame("x_diff" = diff(range(pretty_data[which(pretty_data$ID == runner_ID),]$X)),
          #                             "y_diff" = diff(range(pretty_data[which(pretty_data$ID == runner_ID),]$Y)),
          #                             "ID" = runner_ID)
          #      
          #      ## bind dataframes
          #      x_y_diff = rbind(x_y_diff, runner_df)
          # }
          # 
          # ## extract tracklet IDs
          # ID_to_remove_because_of_less_than_50cm_travel_distance = vector()
          # for (n in 1:nrow(x_y_diff)) {
          #      
          #      ## go through all IDs
          #      runner_ID = x_y_diff$ID[n]
          #      
          #      ## extract row
          #      runner_row = x_y_diff[which(x_y_diff$ID == runner_ID),]
          #      
          #      ## check if tracklet needs removal
          #      if (runner_row$x_diff <= 100 | runner_row$y_diff <= 100) {
          #           
          #           ## list ID in vector
          #           ID_to_remove_because_of_less_than_50cm_travel_distance = c(ID_to_remove_because_of_less_than_50cm_travel_distance,runner_row$ID)  
          #           
          #      }
          #      
          # }
          # 
          # ## remove tracklets
          # for (n in ID_to_remove_because_of_less_than_50cm_travel_distance) pretty_data = pretty_data[-c(which(pretty_data$ID == n)),]
          # 
     ## old marker finding ####
          
          ## identify interesting markers
          
          # ## extract range for all markers and all coordinates
          # range_all_marker_all_coordinates = as.data.frame(apply(raw_data_no_z[,2:ncol(raw_data_no_z)], 2, range, simplify = FALSE))
          # 
          # ## identify all markers which moved at least 3000 units - no idea how meaningful this is (must do for now)
          # moved_marker_col_indices = which(lapply(range_all_marker_all_coordinates,diff) > 3000)
          # 
          # ## extract names from col names from identified markers
          # label_extractor = unlist(regmatches(names(moved_marker_col_indices), gregexpr('\\(?[0-9,.]+', names(moved_marker_col_indices))))
          # label_extracted = as.numeric(gsub('\\(', '-', gsub(',', '', label_extractor)))
          # 
          # ## candidate markers
          # label_extracted
          # 
          # ## identify candidates, where the average step length between two datapoints is realistic (i.e. 10 units)
          # good_candidates = vector()
          # for (n in 1:length(unique(label_extracted))){
          # 
          #      runner_label = unique(label_extracted)[n]
          # 
          #      ## identify coordinates from raw_data
          #      x_col = paste("X",runner_label, sep="")
          #      y_col = paste("Y",runner_label, sep="")
          # 
          #      ## combine coordinates for convenient analysis
          #      test_data = raw_data[,which(names(raw_data) == x_col):which(names(raw_data) == y_col)]
          # 
          #      ## remove 0 values from vector
          #      x_test = test_data[,1]; x_test = x_test[x_test!=0]
          #      y_test = test_data[,2]; y_test = y_test[y_test!=0]
          # 
          #      ## see the differences between points
          #      x_mean_diff = mean(abs(diff(x_test)))
          #      y_mean_diff = mean(abs(diff(y_test)))
          # 
          #      x_length = length(x_test)
          # 
          #      ## store if there are at least 10 datapoints
          #      if (x_length >= 10){
          # 
          #           ## and if the mean differences between data points for at least one axis is is at least 10 units
          #           if (x_mean_diff >= 10 | y_mean_diff >= 10){
          # 
          #                ## store marker
          #                good_candidates[n] = n
          #           }
          #      }
          # }
          # 
          # ## remove NAs from vector
          # good_candidates = good_candidates[which(!is.na(good_candidates))]
          # good_candidates = unique(label_extracted)[good_candidates]
          
          
     ## old stuff below ####  
     ## find tracklet length ####
          
          ## prepare vector
          tracklet_length = vector()
          for (n in 1:length(unique(concise_data$ID))) {
               
               ## run through tracklets
               runner_tracklet = unique(concise_data$ID)[n]
               
               ## get tracklet length
               tracklet_length[n] = length(which(concise_data$ID == runner_tracklet))
               
          }
          concise_data
     ## plot candidate markers ####
          
          ## plot candidates
          for(n in 1:length(unique(good_candidates))) {
               
               ## take marker
               runner_label = unique(good_candidates)[n]
               
               ## identify coordinates from raw_data
               x_col = paste("X",runner_label, sep="")
               y_col = paste("Y",runner_label, sep="")
               
               ## combine coordinates for convenient plotting
               plot_data = raw_data[,which(names(raw_data) == x_col):which(names(raw_data) == y_col)]
               
               ## remove 0 values from vector
               x_plot = plot_data[,1]; x_plot = x_plot[x_plot!=0]
               y_plot = plot_data[,2]; y_plot = y_plot[y_plot!=0]
               
               ## create filename for export
               runner_filename = paste("plot_for_marker_",runner_label,"______time_tracked_",length(x_plot),"_seconds",".png", sep="")
               
               ## create file for export
               png(filename=runner_filename)
               
               ## create plot
               plot(x_plot, y_plot,
                    xlab="X-axis", ylab="Y-axis",
                    xlim=c(-6000,6000), ylim=c(-1000,8000),
                    main=paste("Marker ",runner_label, sep=""),
                    bg="gray50", pch=21)
               
               ## close graphic device
               dev.off()
               
          }
          
     ## use clean data file from Mathias ####
          
          ## load clean data
          mathias_clean_data = read.table(file = paste(raw_data_location, "older stuff/Ant_trails_TerraX.tsv.csv", sep=""), sep = ';', header = FALSE, skip=12)
          
          ## fix number format
          mathias_clean_data_df = as.data.frame(apply(mathias_clean_data, 2, function (x) gsub("\\.","",x)))
          mathias_clean_data_df = as.data.frame(apply(mathias_clean_data_df,2,function (x) paste(substr(x,1,(nchar(x)-3)),".",substr(x,(nchar(x)-2),nchar(x)), sep="")))
          mathias_clean_data_df = apply(mathias_clean_data_df, 2, function(x) as.numeric(x)); mathias_clean_data_df = as.data.frame(mathias_clean_data_df)
          
          ## fix col names
          m=0
          for(n in seq(from=1, to=ncol(mathias_clean_data_df), by=3)) {
               
               ## counter to label the markers
               m = m + 1
               
               ## name each X,y & Z col with the counter m
               names(mathias_clean_data_df)[n:(n+2)] = c(paste("X",m,sep=""), paste("Y",m,sep=""), paste("Z",m,sep=""))
          }
          
          ## export file
          write.table(mathias_clean_data_df, "C:/Users/timte/Downloads/Mathias_data.txt", sep="\t")
          
          ## import file
          clean_df = read.table("C:/Users/timte/Downloads/Mathias_data.txt", sep="\t", header = TRUE)
          
          ## very sloppy way to find drop events in tracklets
          quick_and_dirty_drop_finder = function(data, step_size, min_duration) {
               
               ## make output
               return_list = list()
               
               ## extract step sizes
               x_step = abs(diff(data[,1]))
               y_step = abs(diff(data[,2]))
               
               ## sum step sizes to account for movement on both axis
               sum_step = x_step + y_step
               
               ## identify differences between steps
               instances = abs(round(diff(sum_step),2))
               
               ## make vector with movement / no-movement decision
               instances_movement = vector()
               instances_movement[which(instances <= step_size)] = "no_m"
               instances_movement[which(instances > step_size)] = "m"
               
               ## identify lengths of sequences of non movement (as defined above)
               rle_output = rle(instances_movement)
               
               ## count all events of immobility longer than 3 seconds
               non_movement_event_length = rle_output$lengths[which(rle_output$lengths > min_duration)]
               index_in_rle_output = which(rle_output$lengths %in% non_movement_event_length)
               
               if (length(index_in_rle_output) != 0) {
               
                    ## extract indices
                    index_to_non_movement = vector()
                    for (n in 1:length(index_in_rle_output)) {
                         
                         ## one by one
                         runner = index_in_rle_output[n]
                         
                         ## extract (count) the index from rel_output
                         if (runner != 1) {
                              index_to_non_movement[n] = sum(rle_output$lengths[1:runner-1])
                         } else {
                              index_to_non_movement[n] = 0
                         }
                    }
                    
                    return_list[[1]] = index_to_non_movement
                    return_list[[2]] = non_movement_event_length
                    return_list[[3]] = nrow(data)
                    names(return_list) = c("start_of_drop_events", "lengths_of_drop_events", "length_of_track_in_frames")
               
               } else {
                         
                    return_list[[1]] = NA
                    return_list[[2]] = NA
                    return_list[[3]] = nrow(data)
                    names(return_list) = c("start_of_drop_events", "lengths_of_drop_events", "length_of_track_in_frames")
                    
                    }
               
               return(return_list)
          }
          
          ## process all tracks
          m=0; out_list = list()
          for (n in seq(from=1, to=19, by=3)){
               
               ## to fill list
               m=m+1
               
               ## extract tracks one by one
               runner_df = clean_df[,n:(n+1)]
               runner_df = runner_df[which(runner_df != 0),]
               
               ## find drop events and store in out list
               out_list[[m]] = quick_and_dirty_drop_finder(runner_df, step_size = 3, min_duration = 60)
          }
          names(out_list) = c("Track1", "Track2", "Track3", "Track4", "Track5", "Track6", "Track7")
          
          ## extract lengths
          lengths_of_drop_events_mathias = vector()
          for (n in 1:length(out_list)) lengths_of_drop_events_mathias = c(lengths_of_drop_events_mathias, out_list[[n]]$lengths_of_drop_events)
          
          
          index = out_list[[1]]$start_of_drop_events[2]
          duration = out_list[[1]]$lengths_of_drop_events[2]
          plotter = track_one[index:(index+duration),]; plot(plotter[,1], plotter[,2])
          
          
          
          ## load tracklets
          concise_data = read.table(file = paste(raw_data_location, "older stuff/concise_data.txt", sep=""), sep = "\t", header = TRUE)
          
          ## remove tracklets that are shorter than 1 minute
          rle_output = rle(concise_data$ID)
          subset_concise_data = as.matrix(concise_data[concise_data$ID %in% rle_output$values[which(rle_output$lengths > 600)],])
          
          ## find drop off events in the remaining tracklets
          out_out_list = vector("list", length(unique(subset_concise_data[,"ID"])))
          for (n in 1:length(unique(subset_concise_data[,"ID"]))){
               
               ## one by one
               runner_ID = unique(subset_concise_data[,"ID"])[n]
               
               ## find the drops
               out_out_list[[n]] = quick_and_dirty_drop_finder(subset_concise_data[which(subset_concise_data[,"ID"] == runner_ID),,drop=FALSE],3, 60)
          }
          
          ## extract the lengths of the events
          lengths_of_drop_events = vector()
          for (n in 1:length(out_out_list)) lengths_of_drop_events = c(lengths_of_drop_events, out_out_list[[n]]$lengths_of_drop_events)
          
          ## remove first six (from the L-frame)
          lengths_of_drop_events = lengths_of_drop_events[7:length(lengths_of_drop_events)]
          
          hist(log(lengths_of_drop_events))
          length(lengths_of_drop_events)
          length(which(lengths_of_drop_events >= 600))
          
          
          write.table(lengths_of_drop_events, "C:/Users/timte/Downloads/drop_off_event_lengths.txt", sep="\t")
          
          
          
          ## generate plots
          for (n in 1:(ncol(mattias_clean_data_df)/3)){
               
               ## generate col name
               x_col = paste("X",n,sep="")
               y_col = paste("Y",n,sep="")
               
               ## combine coordinates
               plot_data = mattias_clean_data_df[,which(names(mattias_clean_data_df) == x_col):which(names(mattias_clean_data_df) == y_col)]
               
               ## remove 0 values from vector
               x_plot = as.numeric(plot_data[,1]); x_plot = x_plot[x_plot!=0]
               y_plot = as.numeric(plot_data[,2]); y_plot = y_plot[y_plot!=0]
               
               ## create filename for export
               runner_filename = paste("plot_for_marker_",n,"______time_tracked_",length(x_plot),"_seconds",".png", sep="")
               
               ## create file for export
               png(filename=runner_filename)
               
               ## plot data
               plot(x_plot,y_plot,
                    main=paste("Marker ",n, sep=""),
                    bg="gray50", pch=21)
               
               ## close graphic device
               dev.off()
               
          }
          
     ## console placeholder ####
          
          ## pretty_data time bins
          
          runner = pretty_data[which(pretty_data[,"frame_number"] >= 58000),]
          plotter = runner[which(runner[,"frame_number"] <= 60000),]
          
          plot(plotter[,"X"], plotter[,"Y"],
               xlim=c(-3800,4600), ylim=c(-2200, 8300))
          hist(t(table(stitching_data_first[,"frame_number"])))
          
          size = matrix(data=NA, ncol = 4000, nrow=500)
          object.size(size)
          
          
          
     ## old function, idk ####
          # ## put data in list format
          # stitching_list_first = list(); stitching_list_last = list();
          # for (n in 1:nrow(u_first)) my_list[[n]] = stitching_data_first[n,]; for (n in 1:nrow(u_first)) my_list[[n]] = stitching_data_last[n,]
          # 
          # ## try with different output structure
          # fragment.stitcher.v2 = function (u_first, u_last, u_time_window, u_x_diff, u_y_diff) {
          #      
          #      ## run through stitching_list_last 
          #      
          # }
          