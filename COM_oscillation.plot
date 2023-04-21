# ocilliation of the COM in a stride

# Limb (k) is starts swing at time (t) with parameters (F0), (kv), (Gu) and (Tsw);
# it has velocity (vv), stance time (ttst[k]), previous swing time (ttsw[k]) and duty factor (dutyfactor).

# simulation data saved as: 1 (k),  2 (t),  3 (F0),  4 (kv),  5 (Gu),  6 (Tsw), 7(vv), 8 (ttst[k]), 9 (ttsw[k]), 10 (dutyfactor) 11 (COM x coordinate) 12 (COM y coordinate) 13 (body angle) 14 (number of limbs in swing)

# Set output format and filename
# set terminal pngcairo size 800,600
# set output 'deviation_of_COM.png'

# # Read the data and filter it based on the given conditions
# # plot "<(awk '$23==1 && $15==1' dat)" using 1:($2):($3) with points title 'Key Frames' pointtype 7, \
# #      "dat" using 1:(deviation($2, $3, $2 - (($15 == 1) ? $2 : 0), $3 - (($15 == 1) ? $3 : 0), \
# #      $2 - ((column(-1) == 1) ? $2 : 0), $3 - ((column(-1) == 1) ? $3 : 0))) with lines title 'Deviation of COM'


# # Read the data and filter it based on the given conditions
# key_frames = "<(awk '$23==1 && $15==1' dat)"

# # Define the trajectory function using linear interpolation between key frames
# trajectory(x) = (x >= int_x1) ? (x <= int_x2 ? deviation(x, y, int_x1, int_y1, int_x2, int_y2) : 0) : 0
# set table "key_frames_interp.txt"
# plot key_frames using (int_x1 = $1, int_y1 = $2, int_x2 = $3, int_y2 = $4, x = $1):(y = $2) smooth csplines
# unset table

# # Plot the deviation of all data points from the trajectory
# plot "dat" using 1:(trajectory($2, $3)) with lines title 'Deviation of COM'

reset
set title "Deviation of COM from Trajectory"
set xlabel "Time"
set ylabel "Deviation of COM"

# Function to calculate the distance between a point (x, y) and a line (x1, y1, x2, y2)
distance(x, y, x1, y1, x2, y2) = abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) / sqrt((y2 - y1)**2 + (x2 - x1)**2)

# Function to find keyframes
is_keyframe(x) = (x == 1 ? 1 : 0)

# Process the data
stats "dat" u (is_keyframe(column(23))*is_keyframe(column(15)) == 1 ? $1 : NaN) nooutput

# Find the indices of keyframes in the data file
array keyframe_times[STATS_records]
j = 0
do for [i=1:STATS_records] {
    stats "dat" every ::i-1::i u (is_keyframe(column(23))*is_keyframe(column(15)) == 1 ? $1 : NaN) nooutput
    if (STATS_valid_records > 0) {
        j = j + 1
        keyframe_times[j] = STATS_min
    }
}

# Calculate the deviation
deviation(x, y, t) = distance(x, y, column(2, i), column(3, i), column(2, i+1), column(3, i+1)) for [i=1:STATS_records-1] if (keyframe_times[i] <= t && t <= keyframe_times[i+1])

# Plot the deviation
plot "dat" u 1:(deviation($2, $3, $1)) w lines title "Deviation"



