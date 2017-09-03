# Mouse-tracking Social Simon



This is a complete script for analyzing mouse trajectories collected in the
mouse-tracking Social Simon task. The complete experiment consisted of 4 conditions:

* condition 1: individual Simon task, in which a person had to respond to two colors by moving
right or left and clicking on a response "button"
* condition 2: individual go-nogo version of the Simon task, in which a person had to respond
to only one of the two colors
* condition 3: social Simon task with visual feedback, in which two people carried out the task
together, responding each to one of the colors; the visual feedback consisted in
perceiving the mouse cursor of the co-participant
* condition 4: social Simon task with no visual feedback, in which two people carried out
the task as in the third condition but could only see their own mouse cursor
on the screen

The script goes from data-preprocessing through various ways to analyze mouse-tracking
data that have been applied in the literature. This version of the script
carries out analysis only on the condition 3 as it is of highest interest in this
study. It can be applied directly to analyzing data from condition 4 and with slight
modifications to conditions 1 and 2.


# Data pre-processing

First of all we save the screen dimensions and relevant coordinates used in the experiment.


```r
# save experiment parameter values
screen.width <- 1920
screen.height <- 1080
start.boundary <- 988 # upper boundary of the start button
response.boundary <- 128 # lower boundary of response boxes
stim.boundary <- 820 # the y-coord that needed to be crossed for stimulus to appear
```

We then begin by loading the data and reformatting it for easier processing.


```r
process.joint.mdata <- function(filename) {
    # get all joint data for an experiment and restructure it into long format
    
    mdata <- read.csv(filename)
    names(mdata) <- gsub(x = names(mdata), pattern = "\\.", replacement = "0")
    # coordinates for different people in separate columns
    mdata.long <- reshape(mdata, varying=c(9:dim(mdata)[2]), direction="long", 
                          timevar="order", idvar="trial", sep="_")
    mdata.complete <- mdata.long[complete.cases(mdata.long),]
    # just 2 columns for x and y coordinates
    names(mdata.complete)[11:14] <- c("x_p1", "y_p1", "x_p2", "y_p2")
    molten <- reshape(mdata.complete, varying=c(11:14), timevar="person", 
                      direction="long", sep="_")
    rownames(molten) <- c()
    subid1 <- sort(unique(molten$subid))[1]
    subid2 <- sort(unique(molten$subid))[2]
    
    molten <- mutate(molten, personid=ifelse((person == 'p1'),
                                         subid1, subid2))
    molten$pair <- as.factor(paste(as.character(unique(molten$subid)), 
                                   collapse=""))
    return(molten)
}

# get all experiment data for Experiment 3
path <- paste(getwd(), '/Experiment3', sep="")
file.names <- list.files(path, full.names = TRUE)
exp3 <- do.call(rbind,lapply(file.names, process.joint.mdata))

# rename pair identifiers
pair.levels <- levels(exp3$pair)
exp3$pair <- plyr::mapvalues(exp3$pair, from = pair.levels, 
                             to = c(1:length(pair.levels)))
```

Then we add the coding for independent variables on 

* trial type, which depends on the congruency between the color of the cue and 
its location
* whose turn it was to respond in a given trial, which depends on the color of the cue
* previous trial type and whose turn it was to respond


```r
# add a variable that specifies trial type
exp3 <- mutate(exp3, trial.type = ifelse(cuePos == cueColor, 'congruent', 'incongruent'))
# add a variable that specifies the role
exp3 <- mutate(exp3, role.type=ifelse(((cueColor == 1 & person == 'p1') |
                                        (cueColor == 2 & person == 'p2')),
                                        'active', 'passive' ))

exp3 <- arrange(exp3, pair, personid, trial, order)

# add two variables that specify the trial type and role in the previous trial
get.previous <- function(person.data) {
    needed.cols <- select(person.data, trial, trial.type, role.type, order)
    needed.cols %>% group_by(trial) %>% filter(row_number()==n()) ->
        trial.info.wide
    trial.info.wide$prev.trial = c(NA, trial.info.wide$trial.type[1:(640-1)])
    trial.info.wide$prev.role = c(NA, trial.info.wide$role.type[1:(640-1)])
    
    expanded <- trial.info.wide[rep(1:nrow(trial.info.wide), trial.info.wide$order),]
    
    return(list(expanded$prev.trial, expanded$prev.role))
}

exp3.tb <- data.table(exp3)
exp3 <- as.data.frame(exp3.tb[,  c("prev.trial", "prev.role") := get.previous(.SD), 
                                  by = personid])
rm(exp3.tb)

# select and order the columns
exp3 <- select(exp3, pair, personid, person, block, trial, cuePos, cueColor,
               selectedBox, trial.type, role.type, prev.trial, prev.role,
               startTime, stimTime, order, times, x, y)
```

The y-coordinates are immediately flipped vertically because the package that was
used for collecting the data (Matlab Psychtoolbox) encodes the screen's top left 
as coordinates [0, 0] and therefore y-coordinates grow towards the bottom of the screen
while for ease of analysis we would like them to grow towards the screen's top.

We plot an example trial before any further pre-processing. The complete trajectory 
that is plotted contains all coordinates recorded since the start of the
trial and so also before the participants clicked on the start button located
at the bottom of the screen.


```r
exp3$y <- screen.height - exp3$y
# remap also the relevant boundaries
flipped.start.boundary <- screen.height - start.boundary
flipped.response.boundary <- screen.height - response.boundary
flipped.stim.boundary <- screen.height - stim.boundary

# visualize the data
ggplot(data = filter(exp3, pair==1, trial==20), 
       aes(x=x, y=y, color=person)) +
    geom_path() + xlim(0, screen.width) + ylim(0, screen.height) +
    ggtitle("Complete recorded coordinates of pair 1, trial 20") +
    theme(plot.title = element_text(hjust = 0.5))
```

![](SocialSimon_files/figure-html/flip_y-1.png)<!-- -->

Now we filter out only successful trials, in which participants did not miss any
deadlines and a correct response was given.


```r
# clicked the start button
exp3.complete <- filter(exp3, startTime>0)
# moved as requested
exp3.complete <- filter(exp3.complete, stimTime>0)
# selected the correct response
exp3.complete <- filter(exp3.complete, cueColor==selectedBox)
```

Complete trajectories will be needed for the dynamical part of our analysis at 
the end of this script. For remaining analyses we need to extract
the portions of the trajectories after participants have clicked on the start
button.


```r
# extract trajectory portions after start click
extract.trajectories <- function(person.data, start.border) {
    # print(person.data$trial)
    start.time <- person.data$startTime[1]
    after.startTime <- person.data$times > start.time
    start.index <- min(which(after.startTime==TRUE))
    num.times <- dim(person.data)[1]
    if (person.data$y[start.index] < start.border) {
        # this works only for the person who clicked second
        mdata <- person.data[person.data$times > start.time,]
    } else {
        # for the other person we have to infer...
        before.start <- person.data[1:start.index,]
        own.start.index <- max(which((before.start$y > start.border)==FALSE))+1
        mdata <- person.data[own.start.index:num.times,]
    }
    return(mdata)
}

exp3.complete %>% group_by(pair, personid, trial) %>% 
    do(extract.trajectories(., flipped.start.boundary)) -> 
    exp3.subset

# reset sampling order to start from 0
exp3.subset %>% group_by(personid, trial) %>% 
    do(mutate(., sampling.order=c(1:(max(order) - min(order) + 1)))) -> 
    exp3.subset

# Visualize the data
ggplot(data = filter(exp3.subset, pair==1, trial==20), 
       aes(x=x, y=y, color=person)) +
    geom_path() + xlim(0, screen.width) + ylim(0, screen.height) +
    ggtitle("Coordinates of pair 1, trial 20 after start click") +
    theme(plot.title = element_text(hjust = 0.5))
```

![](SocialSimon_files/figure-html/extract_trajectories-1.png)<!-- -->

For convenience we rescale the coordinates into a standard MouseTracker 
coordinate space, where x is in range [-1, 1] and y in range [0, 1.5].


```r
rescale.space <- function(coordinates, new_min, new_max, old_min, old_max){
    # Rescales space to given x- and y-limits and centers on the origin.
    # :param coordinates: (array) array of coordinates to rescale
    # :param new_min: (int) new minimum limit
    # :param new_max: (int) new maximum limit
    # :param old_range: (int) original range of the data (screen width or height)
    # :return: space-normalized array
    new_range <- new_max - new_min
    old_range <- old_max - old_min
    rescaled <- ((((coordinates - old_min) * new_range) / old_range) + new_min)
    return(rescaled)
}

exp3.subset$x.scaled <- rescale.space(exp3.subset$x, -1, 1, 0, screen.width)
exp3.subset$y.scaled <- rescale.space(exp3.subset$y, 0, 1.5, 0, screen.height)

scaled.start.boundary <- rescale.space(flipped.start.boundary, 0, 1.5, 0, 
                                       screen.height)
scaled.response.boundary <- rescale.space(flipped.response.boundary, 0, 1.5, 0, 
                                          screen.height)
scaled.stim.boundary <- rescale.space(flipped.stim.boundary, 0, 1.5, 0, 
                                      screen.height)

# Visualize the data
ggplot(data = filter(exp3.subset, pair==1, trial==20), 
       aes(x=x.scaled, y=y.scaled, color=person)) +
    geom_path() + xlim(-1, 1) + ylim(0, 1.5) +
    ggtitle("Space rescaled coordinates of pair 1, trial 20") +
    theme(plot.title = element_text(hjust = 0.5))
```

![](SocialSimon_files/figure-html/space_rescaling-1.png)<!-- -->

Now we align all the trajectories to the common [0, 0] origin and timestamps
to start at 0.


```r
normalize.space <- function(coordinates) {
    new.coord <- coordinates - coordinates[1]
    return(new.coord)
}

get.stim.time <- function(stim.time, times) {
    return(stim.time - times[1])    
}

exp3.subset %>% group_by(pair, trial, person) %>% 
    do(mutate(., x.aligned = normalize.space(x.scaled), 
              y.aligned = normalize.space(y.scaled),
              t.aligned = normalize.space(times),
              stim.aligned = get.stim.time(stimTime, times))) -> 
    exp3.aligned

#given that we have no data regarding previous trial for trial 1, remove them
exp3.aligned <- exp3.aligned[complete.cases(exp3.aligned),]

# Visualize data
ggplot(data = filter(exp3.aligned, pair==1, trial==20), 
       aes(x=x.aligned, y=y.aligned, color=person)) +
    geom_path() + xlim(-1, 1) + ylim(0, 1.5) +
    ggtitle("Space aligned coordinates of pair 1, trial 20") +
    theme(plot.title = element_text(hjust = 0.5))
```

![](SocialSimon_files/figure-html/align_trajectories-1.png)<!-- -->

At this point we can already visualize all trajectories of all pairs.


```r
# Visualize all pairs
ggplot(data = exp3.aligned, aes(x=x.aligned, y=y.aligned, 
                                color=person, group=trial)) +
    geom_path() + xlim(-1, 1) + ylim(0, 1.5) + facet_wrap(~pair) +
    ggtitle("Trajectories for all pairs for all trials") +
    theme(plot.title = element_text(hjust = 0.5))
```

![](SocialSimon_files/figure-html/all_trajectories-1.png)<!-- -->

We see that all but one pair seem to divide the screen space between each other
by moving mostly directly towards their assigned response box and avoiding the center.
However, one particular pair (number 7) in condition 3 (with visual feedback) has 
mostly upward moving trajectories. We can further explore here
whether the joint upward motion is induced by one of the participants or happens
immediately on both sides.


```r
# plot the first 5 trials of aligned trajectories
p7.early <- filter(exp3.aligned, pair==7, trial %in% c(2:7))
ggplot(data = p7.early, 
       aes(x=x.aligned, y=y.aligned, group=trial, color=role.type)) +
    geom_path() + xlim(-1, 1) + ylim(0, 1.5) + 
    facet_grid(trial~person)
```

![](SocialSimon_files/figure-html/aberrant_pair-1.png)<!-- -->

The plot presents trajectories in the first couple of trials (2, 3, 5, 6, 7) of two
participants from pair 7. The trajectories are colored red for when participant's
role was "active", i.e. the cue that appeared had their assigned color and blue for
when their role was "passive", i.e. their task in such a trial was to not press
the response button and simply go back to the start position.

From the plots it would appear that one of the members of this aberrant couple
(person 2) adopts a "move upward" strategy from the start, independently of whether
it is their turn to respond or not. The other person in that couple seems to copy
the co-actor's movements on passive trials and make large discrete errors on the trials
in which it is actually their turn to respond. Why this particular couple behaves
in this manner is unfortunately unknown and might indicate some individual differences,
the feeling of jointness experienced by this couple or conscious strategies that 
people adopt in such a task.

As the next pre-processing step we will flip all trajectories to one side.
This ensures that every trajectory starts at the bottom of the coordinate system 
and ends in the top right corner. It is done to obtain comparable trajectory 
measures.


```r
flip.x <- function(coordinates) {
    n.coord <- length(coordinates)
    n.neg <- sum(coordinates<0)
    
    if (coordinates[n.coord] < -0.1 || n.neg > (0.5*n.coord)) {
        coordinates <- coordinates * -1
    }
    return(coordinates)
}

exp3.aligned %>% group_by(pair, trial, person) %>% 
    do(mutate(., x.flipped = flip.x(x.aligned))) -> 
    exp3.aligned

# also convert times into ms
exp3.aligned$t.aligned <- exp3.aligned$t.aligned*1000

# remove irrelevant columns
exp3.aligned %>% 
    select(-c(startTime, stimTime, order, times, x, y, 
              x.scaled, y.scaled, x.aligned)) -> 
    exp3.aligned.clean
```

As a last step, a look at sampling rate distribution to see whether it reveals 
any outliers that could indicate missing data or wrong recording.


```r
exp3.aligned.clean %>% group_by(personid, trial) %>% 
    do(as.data.frame(diff(.$t.aligned, lag=1))) -> 
    tdiff
colnames(tdiff)[3] <- 'intervals'
rate.mean <- mean(tdiff$intervals)
num.outliers <- sum(tdiff$intervals > 3 * sd(tdiff$intervals) + rate.mean)
mean.outliers <- mean(tdiff$intervals[which(tdiff$intervals > 3 * sd(tdiff$intervals) + 
                          mean(tdiff$intervals))])

# y coordinates below -0.5 are a result of wrong recording of the start press
# sum(exp3.aligned.clean$y.aligned < -0.1)
exp3.aligned.clean %>% ungroup() %>%
    filter(y.aligned < -0.1) %>%
    select(pair, trial) %>% unique() ->
    yglitch
yglitch <- as.data.frame(yglitch)
yglitch$pair <- as.integer(as.character(yglitch$pair))

for (i in 1:dim(yglitch)[1]) {
    exp3.aligned.clean <- filter(exp3.aligned.clean, 
                                 !(pair == yglitch[i,1] & trial == yglitch[i,2]))
}
```
We see that the mean of the overall sampling rate is 10.87 ms, which is to 
be expected given that the experiment sampling rate was set to 92 Hz. However, in 
41 instances the time interval between one sampling point and the next is higher 
than 3 standard deviations away from the mean rate. The mean of all those points 
is around 207.71 ms and indicates an interruption in data recording.
In order to facilitate binned data analysis, we remove trials that contain these
large large deviations.


```r
# identify which person the data comes from
pmissing <- unique(tdiff$personid[which(tdiff$intervals > 2*sd(tdiff$intervals) +
                                    mean(tdiff$intervals))])
# which trial
tmissing <- unique(tdiff$trial[which(tdiff$intervals > 2*sd(tdiff$intervals) +
                                     mean(tdiff$intervals))])

exp3.aligned.clean <- subset(exp3.aligned.clean, 
                    !(personid %in% pmissing[1:2] & trial %in% tmissing[1:10]))
exp3.aligned.clean <- subset(exp3.aligned.clean, 
                            !(personid %in% pmissing[3:4] & trial %in% tmissing[11:21]))
```

With the sampling timing issues fixed, we can examine reaction time outliers and
add another variable to the data, which indicates whether the trial was fast or slow.


```r
get.RT <- function(time_coordinates) {
    n.coord <- length(time_coordinates)
    return(time_coordinates[n.coord])
}

exp3.aligned.clean %>% group_by(pair, trial, person) %>% 
    do(mutate(., total.time = get.RT(t.aligned))) %>%
    ungroup() ->
    exp3.aligned.clean

exp3.aligned.clean %>% group_by(personid) %>%
    mutate(., median.time = median(total.time)) %>%
    ungroup() ->
    exp3.aligned.clean
    
ggplot(exp3.aligned.clean, aes(x=total.time)) + geom_histogram(binwidth=100) +
    ggtitle("Histogram of all trial durations")
```

![](SocialSimon_files/figure-html/trial_speed-1.png)<!-- -->

```r
# remove particularly slow trials
outlierRT <- mean(exp3.aligned.clean$total.time) + 
    3*sd(exp3.aligned.clean$total.time)
exp3.aligned.clean %>% group_by(pair, trial, person) %>% 
    filter(total.time < outlierRT) -> 
    exp3.aligned.clean

ggplot(exp3.aligned.clean, aes(x=total.time)) + geom_histogram(binwidth=100) +
    ggtitle("Histogram of trial durations with particularly long trials removed")
```

![](SocialSimon_files/figure-html/trial_speed-2.png)<!-- -->

```r
# add 
midRT <- median(exp3.aligned.clean$total.time)
exp3.aligned.clean %>% group_by(pair, trial, person) %>% 
    do(mutate(., trial.speed = ifelse(total.time <= midRT,
                                      'fast', 'slow'))) -> 
    exp3.aligned.clean
```

For the analysis we will use a mousetrap package. Accordingly, the next step is 
to transform the data into a mousetrap object.


```r
mtdata <- mt_import_long(exp3.aligned.clean, 
                        xpos_label = 'x.flipped', 
                        ypos_label = 'y.aligned', 
                        timestamps_label = 't.aligned', 
                        mt_id_label = c('personid','trial'), 
                        mt_seq_label = 'sampling.order', 
                        reset_timestamps = FALSE)
```

As a final pre-processing step, we perform time normalization on the data, in 
which the times and coordinates are linearly interpolated so that each trajectory 
contains the same number of recorded points, typically this is set to 101 points.


```r
mtdata <- mt_time_normalize(mtdata)
```


# Dependent variables calculation

There is a number of measures that can be calculated on the basis of raw time
and normalized trajectories. 
First, we retrieve trajectory derivatives (velocity, acceleration) and angles 
from the raw time data. Next, we calculate a variety of measures on normalized
data. 


```r
# get mt measures
mtdata <- mt_derivatives(mtdata, use="trajectories")
mtdata <- mt_angles(mtdata, use="trajectories")
mtdata <- mt_measures(mtdata, use="trajectories", save_as = "measures")
# TODO implement distance calculation from response and foil coordinates
# TODO we can calculate
# max_angle:      max angle towards incorrect response
# init_angle:     initial angle after the start of the movement
# angle_flips:    the amount of change in angles

# this returns
# a) vertical-based angles (angle_v), where positive values indicate a movement 
# to the left of the vertical, negative - to the right 
# b) point-based angles (angle_p) that indicate changes of movement within three
# consecutive time steps; a value of pi (= 3.14...) (for radians) or 180 (for 
# degrees) indicates a constant movement direction, a value of 0 a complete reversal
# Freeman et al suggest calculating angle relative to the start button.
```

These are the measures that are returned:

| Measure | Definition |
|:--------|:-----------------------------------------------|
|xpos_max|   Maximum x-position|
|xpos_min|   Minimum x-position|
|ypos_max|   Maximum y-position|
|ypos_min|   Minimum y-position|
|MAD|        Signed Maximum absolute deviation from the direct path connecting start and end point of the trajectory (straight line). If the MAD occurs above the direct path, this is denoted by a positive value; if it occurs below, by a negative value.|
| MAD_time|   Time at which the maximum absolute deviation was reached first|
| MD_above|   Maximum deviation above the direct path|
| MD_above_time|  Time at which the maximum deviation above was reached first|
| MD_below|   Maximum deviation below the direct path|
| MD_below_time|  Time at which the maximum deviation below was reached first|
| AD|         Average deviation from direct path|
| AUC|        Area under curve, the geometric area between the actual trajectory and the direct path where areas below the direct path have been subtracted|
| xpos_flips|         Number of directional changes along x-axis (exceeding the distance specified in flip_threshold)|
| ypos_flips|         Number of directional changes along y-axis (exceeding the distance specified in flip_threshold)|
| xpos_reversals|     Number of crossings of the y-axis|
| ypos_reversals|     Number of crossings of the x-axis|
| RT|                 Response time, time at which tracking stopped|
| initiation_time|    Time at which first mouse movement was initiated|
| idle_time|    Total time without mouse movement across the entirety of the trial|
| hover_time|         Total time of all periods without movement in a trial (whose duration exceeds the value specified in hover_threshold)|
| hovers|             Number of periods without movement in a trial |
| total_dist|     Total euclidean distance covered by the trajectory|
| vel_max|       Maximum velocity|
| vel_max_time|   Time at which maximum velocity occurred first|
| vel_min|   Minimum velocity|
| vel_min_time|   Time at which minimum velocity occurred first|
| acc_max|        Maximum acceleration|
| acc_max_time|   Time at which maximum acceleration occurred first|
| acc_min|        Minimum acceleration|
| acc_min_time|   Time at which minimum acceleration occurred first|

Out of these measures, the most widely used are measures of the curvature of mouse
trajectories, i.e. MAD, AD and AUC. A variety of time-based measures allows us
to examine movement duration in different stages of the trajectory. Finally,
coordinate flips and reversals are a proxy for movement complexity (which we
later examine also with entropy-based measure).


Our data contains trials in which a participant's role was 'active', i.e.,
it was their turn to respond to the cue and trials in which the role was 'passive',
i.e., their task was to refrain from responding. It is reasonable to assume
that different cognitive processes are at play in these types of trials and that
different trajectories should result. Therefore, for further analysis we split 
the data into two groups: active and passive data.


```r
activedata <- mt_subset(mtdata, role.type=='active')
passivedata <- mt_subset(mtdata, role.type=='passive')
```


# Exploratory Analysis

First we should visually examine trajectories averaged across trials
for each participant and across participants. They are averaged separately
for each independent variable of interest, which would allow us to reveal
different patterns depending on the variables we decide to include (if there are
such patterns in the data obviously).

Our independent variables that we will examine are:

* trial type: congruent or incongruent
* previous trial type: whether trial at time t-1 was congruent or not
* previous role type: whether in the previous trial participant had to act


```r
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
# plot average trajectory for active trials, by trial type
p1 <- mt_plot_aggregate(activedata, use = "tn_trajectories", points=TRUE, 
                  x = "xpos", y = "ypos", color = "trial.type", 
                  subject_id = "personid") +
    theme(legend.position=c(.8,.2)) +
    xlab("x coordinate (px)") + ylab("y coordinate (px)") +
    ggtitle("Average trajectories for active trials by trial type")

# plot average trajectory for passive trials, by trial type
p2 <- mt_plot_aggregate(passivedata, use = "tn_trajectories", points=TRUE, 
                  x = "xpos", y = "ypos", color = "trial.type", 
                  subject_id = "personid") +
    theme(legend.position=c(.8,.2)) +
    xlab("x coordinate (px)") + ylab("y coordinate (px)") +
    ggtitle("Average trajectories for passive trials by trial type")

# trajectories by previous trial type
p3 <- mt_plot_aggregate(activedata, use = "tn_trajectories", points=TRUE, 
                  x = "xpos", y = "ypos", color = "prev.trial", 
                  subject_id = "personid") +
    theme(legend.position=c(.8,.2))+
    xlab("x coordinate (px)") + ylab("y coordinate (px)") +
    ggtitle("Average trajectories for active trials by previous trial type")
p4 <- mt_plot_aggregate(passivedata, use = "tn_trajectories", points=TRUE, 
                  x = "xpos", y = "ypos", color = "prev.trial", 
                  subject_id = "personid") +
    theme(legend.position=c(.8,.2))+
    xlab("x coordinate (px)") + ylab("y coordinate (px)")+
    ggtitle("Average trajectories for passive trials by previous trial type")

# trajectories by previous role
p5 <- mt_plot_aggregate(activedata, use = "tn_trajectories", points=TRUE, 
                  x = "xpos", y = "ypos", color = "prev.role", 
                  subject_id = "personid") +
    theme(legend.position=c(.8,.2))+
    xlab("x coordinate (px)") + ylab("y coordinate (px)")+
    ggtitle("Average trajectories for active trials by previous role")
p6 <- mt_plot_aggregate(passivedata, use = "tn_trajectories", points=TRUE, 
                  x = "xpos", y = "ypos", color = "prev.role", 
                  subject_id = "personid") +
    theme(legend.position=c(.8,.2))+
    xlab("x coordinate (px)") + ylab("y coordinate (px)")+
    ggtitle("Average trajectories for passive trials by previous role")

multiplot(p1, p2, p3, p4, p5, p6, cols=2)
```

![](SocialSimon_files/figure-html/explo_plots-1.png)<!-- -->

Already from these plots it would seem that neither of the independent variables
affects the shape of movement trajectory. However, we will persevere. Instead of 
plotting both x,y-coordinates, we could focus our attention on just
the x plane, since it is the one of bigger relevance to the task.


```r
# Plot aggregate trajectories with standard errors
# (note that mean_se does not take into account within subjects design)
avg.tn.trajectories <- mt_aggregate_per_subject(activedata,
  use="tn_trajectories", use2_variables="trial.type",subject_id="personid")
p1 <- ggplot(avg.tn.trajectories, aes(x=steps,y=xpos,group=trial.type))+
    stat_summary(geom = "ribbon",fun.data=mean_se,alpha=.2)+
    geom_line(aes(color=trial.type),stat="summary",fun.y="mean")+
    scale_color_brewer(type="qual",palette = "Set1" )+
    theme(legend.position=c(.2,.8), plot.title = element_text(size = 10)) +
    ggtitle("Average x dimension by trial type")
    
avg.tn.trajectories <- mt_aggregate_per_subject(activedata,
  use="tn_trajectories", use2_variables="prev.trial",subject_id="personid")
p2 <- ggplot(avg.tn.trajectories, aes(x=steps,y=xpos,group=prev.trial))+
    stat_summary(geom = "ribbon",fun.data=mean_se,alpha=.2)+
    geom_line(aes(color=prev.trial),stat="summary",fun.y="mean")+
    scale_color_brewer(type="qual",palette = "Set1" )+
    theme(legend.position=c(.2,.8), plot.title = element_text(size = 10)) +
    ggtitle("Average x dimension by previous trial")

avg.tn.trajectories <- mt_aggregate_per_subject(activedata,
  use="tn_trajectories", use2_variables="prev.role",subject_id="personid")
p3 <- ggplot(avg.tn.trajectories, aes(x=steps,y=xpos,group=prev.role))+
    stat_summary(geom = "ribbon",fun.data=mean_se,alpha=.2)+
    geom_line(aes(color=prev.role),stat="summary",fun.y="mean")+
    scale_color_brewer(type="qual",palette = "Set1" )+
    theme(legend.position=c(.2,.8), plot.title = element_text(size = 10))+
    ggtitle("Average x dimension by previous role")

multiplot(p1, p2, p3, cols=3)
```

![](SocialSimon_files/figure-html/x_plots-1.png)<!-- -->

Rather than looking at grand averages, we can also look at averages within pairs.
For simplicity we only look at trajectories by current trial type.


```r
pairs <- sort(rep(c(1:11), 2*2*101))
avg.tn.trajectories.active <- mt_aggregate_per_subject(activedata,
  use="tn_trajectories", use2_variables="trial.type", subject_id="personid")
avg.tn.trajectories.active$couple <- as.factor(pairs)

avg.tn.trajectories.active <- mutate(avg.tn.trajectories.active, 
                              person=ifelse(personid %% 2 != 0,
                                        'p1', 'p2' ))

ggplot(avg.tn.trajectories.active, aes(x=xpos, y=ypos, 
                     linetype=trial.type, color=person)) + geom_path() + 
    facet_wrap(~couple)
```

![](SocialSimon_files/figure-html/avg_pairs-1.png)<!-- -->

Also here, even if there is variability in trajectories between people, it does 
not seem like congruency of the trial within each person leads to different 
movement path.

We can do the same for the passive data.


```r
avg.tn.trajectories.passive <- mt_aggregate_per_subject(passivedata,
  use="tn_trajectories", use2_variables="trial.type", subject_id="personid")
avg.tn.trajectories.passive$couple <- as.factor(pairs)

avg.tn.trajectories.passive <- mutate(avg.tn.trajectories.passive, 
                              person=ifelse(personid %% 2 != 0,
                                        'p1', 'p2' ))

ggplot(avg.tn.trajectories.passive, aes(x=xpos, y=ypos, 
                     linetype=trial.type, color=person)) + geom_path() + 
    facet_wrap(~couple)
```

![](SocialSimon_files/figure-html/avg_pairs_passive-1.png)<!-- -->

From the examination of passive trials it seems that people tend to complete
the movement, even if they could return to the starting position as soon as
it becomes clear it is not their turn to respond.
We can confirm this by also plotting a histogram of maximum y coordinate
reached in passive trials. The thresholds indicated on the plot are the locations
of the start boundary, the y-coordinated that had to be crossed in order for the 
cue to appear and the lower response box boundary. As can be seen from the plot,
the majority of trajectories goes beyond that last threshold.


```r
cuts <- data.frame(Thresholds=c('start', 'stimulus', 'response'), 
                   vals=c(scaled.start.boundary, scaled.stim.boundary, 
                          scaled.response.boundary))
# we need to use data before alignment was performed
exp3.aligned %>% group_by(pair, trial, person) %>%
    filter(row_number()==n()) %>% select(y.scaled) ->
        max.ys
ggplot(max.ys, aes(x=y.scaled)) + 
    geom_histogram(binwidth=.05) + 
    geom_vline(data=cuts, aes(xintercept=vals, color=Thresholds))
```

![](SocialSimon_files/figure-html/passive_hist-1.png)<!-- -->

In addition to plotting time-normalized trajectories in full, we can bin them into
several intervals.


```r
# average time-normalized trajectories into 4 bins
activedata <- mt_average(activedata, use="tn_trajectories", 
                     save_as="av_tn_trajectories",
                     av_dimension="steps",
                     intervals = c(0.5,33.5,67.5,101.5))

# to view the data
# mt_aggregate(activedata,use = "av_tn_trajectories",
#   use2_variables = "trial.type", subject_id = "personid")

mt_plot_aggregate(activedata, use = "av_tn_trajectories", points=TRUE,
                  color = "trial.type", subject_id = "personid") +
    ggtitle("Active trajectory coordinates in binned normalized time")
```

![](SocialSimon_files/figure-html/binned_tnplots-1.png)<!-- -->

Some analyses require retaining trajectories in raw time. In this case, we
decide how many raw time bins to create between 0 ms and some cutoff (e.g., 
1500 ms) and then create a number of raw time steps. Thus, each step 
(i.e., coordinate pair) of a trajectory reflects the location of the mouse during 
some raw time bin (e.g., 500-600 ms if bins are 100 ms wide).

Such raw binned data can also be plotted although it seems to work better if 
data for slow and fast trials is separated.


```r
# split data in slow and fast trials
slow.activedata <- mt_subset(activedata, trial.speed=='slow')
fast.activedata <- mt_subset(activedata, trial.speed=='fast')

slow.activedata <- mt_average(slow.activedata,
                     use="trajectories", save_as = "av_trajectories",
                     interval_size = 300, max_interval = 900)
mt_plot_aggregate(slow.activedata, use = "av_trajectories", points=TRUE,
                  color = "trial.type", subject_id = "personid") +
    theme(legend.position=c(.8,.2)) +
    xlab("x coordinate (px)") + ylab("y coordinate (px)") +
    ggtitle("Slow active trajectory coordinates in binned raw time")
```

![](SocialSimon_files/figure-html/binned_plots-1.png)<!-- -->

```r
fast.activedata <- mt_average(fast.activedata,
                       use="trajectories", save_as = "av_trajectories",
                       interval_size = 300, max_interval = 1800)
mt_plot_aggregate(fast.activedata, use = "av_trajectories", points=TRUE, 
                  color = "trial.type", subject_id = "personid") +
    theme(legend.position=c(.8,.2)) +
    xlab("x coordinate (px)") + ylab("y coordinate (px)") +
    ggtitle("Fast active trajectory coordinates in binned raw time")
```

![](SocialSimon_files/figure-html/binned_plots-2.png)<!-- -->

Similar analysis can be performed not on coordinates but on movement velocity.

The literature suggests that stronger competition between response options should 
be characterized by an initial decreased velocity as competing choices inhibit 
each other, followed by an increase in velocity once the system converges upon
a decision and the inhibition is alleviated. Thus, analyzing velocity data can 
allow for inferences about when commitments to a particular response are made.

Here we plot velocity profiles in binned raw time, together with approximate
time in which cue appeared.


```r
scaled.stim.times <- mtdata$data$stim.aligned * 1000
stim.mean <- mean(scaled.stim.times)
stim.sd <- sd(scaled.stim.times)
stim.lower <- stim.mean - stim.sd
stim.higher <- stim.mean + stim.sd

cuts.velocity <- data.frame(Thresholds=c('stim.lower', 'stim.mean', 'stim.higher'), 
                   vals=c(stim.lower, stim.mean, 
                          stim.higher))
cuts.velocity$stat <- c('sd', 'mean', 'sd')
cuts.velocity
```

```
##    Thresholds     vals stat
## 1  stim.lower 181.7183   sd
## 2   stim.mean 301.6412 mean
## 3 stim.higher 421.5640   sd
```

```r
# velocity in active data
activedata <- mt_average(activedata,
                     use="trajectories", save_as = "av_trajectories",
                     interval_size = 200, max_interval = 1800)

mt_plot_aggregate(activedata, use = "av_trajectories", 
                  points=TRUE, x = "timestamps", y = "vel", 
                  color = "trial.type", subject_id = "personid") +
    theme(legend.position=c(.75,.85), legend.box='horizontal')+
    xlab("time") + ylab("velocity") +
    geom_vline(data=cuts.velocity, aes(xintercept=vals, linetype=stat)) +
    scale_linetype_manual(name="Stimulus Time Stats", 
                          values=c("solid", "longdash")) +
    ggtitle("Active velocity profiles in binned raw time")
```

![](SocialSimon_files/figure-html/active_vel-1.png)<!-- -->

```r
# velocity in passive data
passivedata <- mt_average(passivedata,
                       use="trajectories", save_as = "av_trajectories",
                       interval_size = 200, max_interval = 1800)
mt_plot_aggregate(passivedata, use = "av_trajectories", 
                  points=TRUE, x = "timestamps", y = "vel", 
                  color = "trial.type", subject_id = "personid") +
    theme(legend.position=c(.75,.85), legend.box='horizontal')+
    xlab("time") + ylab("velocity") +
    geom_vline(data=cuts.velocity, aes(xintercept=vals, linetype=stat)) +
    scale_linetype_manual(name="Stimulus Time Stats", 
                          values=c("solid", "longdash")) +
    ggtitle("Passive velocity profiles in binned raw time")
```

![](SocialSimon_files/figure-html/active_vel-2.png)<!-- -->


# Statistical analysis

Exploratory analysis suggests that apart from an aberrant couple 7, there do not
seem to be any differences in trajectories or velocity profiles either in active
or passive data, nor when we take into account independent variables of interest
(current trial type, previous trial type, previous role).

We will, however, perform statistical analyses to try to confirm these observations.

Given the negative conclusions from the exploratory analysis, we will focus on 
the difference between trajectories depending on current trial type (congruent vs
incongruent) and ignore the effect of previous trial for the time being.


## Testing trajectories directly

One way to test whether trajectories differ in different conditions is by
examining directly coordinates of interest. Since the x-coordinate plane is
typically thought to be more relevant to the Simon task, we focus on analyzing
x coordinates.


### Paired t-tests on coordinates

One  approach is to use 101 paired-samples t tests to compare the  x-coordinate 
of participants' mean trajectories for two conditions at each individual time step. 


```r
xpos.ttests <- 
    with(avg.tn.trajectories.active,
         sapply(unique(steps),function(i){
             t.test(xpos[trial.type=="congruent" & steps==i],
                    xpos[trial.type=="incongruent" & steps==i],
                    paired = TRUE)$p.value})
    )

# Retrieve all significant t-tests
which.sig <- which(xpos.ttests<.05)
# Number of significant t-tests
n.sig <- sum(xpos.ttests<.05,na.rm=TRUE)
# Number of adjacent significant t-tests (minus 1)
# table(diff(which(xpos_t_tests<.05)))
```

The test revealed a sequence of 0 significant t-tests on the difference
between x-coordinates in congruent and incongruent trials. In order to determine
what is the minimum number of significant t-tests that qualifies as a pattern, a
bootstrapping procedure would be required. However, in previous research 8 was
the minimum so this number seems rather low by comparison.

## Anovas on binned trajectories

Other than looking at particular coordinates, we can also run tests on binned
trajectories (both normalized and raw time), that we have also plotted above.
In this case the analysis we perform is repeated measures 3 (bins) by 2 (trial
type) ANOVA. 


```r
avg.tn.trajectory.bins <-  mt_aggregate_per_subject(activedata, 
                                                    use="av_tn_trajectories",
                                                    use2_variables="trial.type",
                                                    subject_id="personid")
avg.tn.trajectory.bins$bin <- factor(avg.tn.trajectory.bins$steps)

aov_ez(data=avg.tn.trajectory.bins,
       id="personid", dv="xpos", within = c("trial.type", "bin"),
       anova_table = list(es=c("ges", "pes"), correction=c("GG")))
```

```
## Anova Table (Type 3 tests)
## 
## Response: xpos
##           Effect          df  MSE          F    ges pes p.value
## 1     trial.type       1, 21 0.00       0.24 <.0001 .01     .63
## 2            bin 1.39, 29.17 0.02 355.38 ***    .80 .94  <.0001
## 3 trial.type:bin 1.64, 34.51 0.00       2.53  .0003 .11     .10
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '+' 0.1 ' ' 1
## 
## Sphericity correction method: GG
```

The result tells us that there is an obvious significant difference in x positions
in different time bins (expected given the nature of the task). However, there
is not significant difference by trial type.


```r
avg.slow.trajectory.bins <- mt_aggregate_per_subject(slow.activedata,
                                                     use="av_trajectories",
                                                     use2_variables="trial.type",
                                                     subject_id="personid")
avg.slow.trajectory.bins$bin <- factor(avg.slow.trajectory.bins$timestamps)

avg.fast.trajectory.bins <- mt_aggregate_per_subject(fast.activedata,
                                                     use="av_trajectories",
                                                     use2_variables="trial.type",
                                                     subject_id="personid")
avg.fast.trajectory.bins$bin <- factor(avg.fast.trajectory.bins$timestamps)

# check that each subject has observations in each category
# table(table(avg.slow.trajectory.bins$personid))
# table(table(avg.fast.trajectory.bins$personid))

aov_ez(data=avg.slow.trajectory.bins, id="personid", dv="xpos", 
       within = c("trial.type","bin"),
       anova_table = list(es=c("ges","pes"), correction=c("GG")))
```

```
## Anova Table (Type 3 tests)
## 
## Response: xpos
##           Effect          df  MSE          F   ges pes p.value
## 1     trial.type       1, 21 0.00       0.90 .0002 .04     .35
## 2            bin 1.64, 34.38 0.02 212.87 ***   .63 .91  <.0001
## 3 trial.type:bin 1.49, 31.28 0.00       2.38  .001 .10     .12
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '+' 0.1 ' ' 1
## 
## Sphericity correction method: GG
```

```r
aov_ez(data=avg.fast.trajectory.bins, id="personid", dv="xpos", 
       within = c("trial.type","bin"),
       anova_table = list(es=c("ges","pes"), correction=c("GG")))
```

```
## Anova Table (Type 3 tests)
## 
## Response: xpos
##           Effect          df  MSE          F    ges pes p.value
## 1     trial.type       1, 20 0.00       0.61 <.0001 .03     .45
## 2            bin 1.41, 28.12 0.02 323.64 ***    .84 .94  <.0001
## 3 trial.type:bin 1.85, 37.05 0.00       1.47  .0003 .07     .24
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '+' 0.1 ' ' 1
## 
## Sphericity correction method: GG
```

```r
# TODO: also velocity
# also on eulidean distance
```

Similarly in raw binned data we find no evidence for significant effect of trial
type.


## Anovas on dependent measures

Having dealt with coordinates, we can move on to looking at dependent measures
typically examined in movement trajectories. That is, the question here is whether
any of the particular measures that summarize trajectories (e.g. maximum x or y
positions, area under curve, x flips) differs between congruent and incogruent
trials, separately for active and passive data.


```r
aggregated.measures <- mt_aggregate_per_subject(mtdata, 
                                                use = 'measures',
                                                subject_id = "personid", 
                                                use2_variables = c("trial.type", 
                                                                   "role.type"))

aggregated.measures %>%
    group_by(role.type, trial.type) %>%
    select(-personid) %>%
    summarize_all(.funs = c("mean","sd")) ->
    measures.summary

# paired t-tests on measures separately for active and passive data
measures <- colnames(aggregated.measures)[4:length(colnames(aggregated.measures))]
n.measures <- length(measures)
outcomes <- matrix(nrow=n.measures*2, ncol=9)

for (i in seq(1, n.measures)) {
    m <- measures[i]
    active <- t.test(eval(as.symbol(m)) ~ trial.type, 
                     data=filter(aggregated.measures, role.type=='active'), 
                     paired=TRUE)
    r <- round(sqrt(active$statistic^2 / (active$statistic^2 + active$parameter)), 3)
    new.row <- c(m, 'active', round(active$statistic, 3), active$parameter, 
                 round(active$p.value, 5), round(active$conf.int[1], 3), 
                 round(active$conf.int[2], 3), round(active$estimate, 3), r)
    outcomes[2*i-1,] <- new.row

    passive <- t.test(eval(as.symbol(m)) ~ trial.type, 
                      data=filter(aggregated.measures, role.type=='passive'), 
                      paired=TRUE)
    r <- round(sqrt(active$statistic^2 / (active$statistic^2 + active$parameter)), 3)
    new.row <- c(m, 'passive', round(passive$statistic, 3), passive$parameter, 
             round(passive$p.value, 5), round(passive$conf.int[1], 3), 
             round(passive$conf.int[2], 3), round(passive$estimate, 3), r)
    outcomes[2*i,] <- new.row
}

outcomes <- data.frame(outcomes)
colnames(outcomes) <- c('measure', 'type', 't', 'df', 'p', 'lower.conf', 
                        'upper.conf', 'estimate', 'effect')

for(i in c(3:ncol(outcomes))) {
    outcomes[,i] <- as.numeric(as.character(outcomes[,i]))
}


significant <- filter(outcomes, p<0.05)
kable(significant)
```



measure           type            t   df         p   lower.conf   upper.conf   estimate   effect
----------------  --------  -------  ---  --------  -----------  -----------  ---------  -------
ypos_max          active     -2.141   21   0.04416       -0.006        0.000     -0.003    0.423
MD_above          passive    -2.260   21   0.03459       -0.015       -0.001     -0.008    0.289
MD_above_time     active     -2.929   21   0.00802      -33.459       -5.675    -19.567    0.539
ypos_flips        passive    -2.218   21   0.03773       -0.091       -0.003     -0.047    0.280
initiation_time   passive     3.320   21   0.00326        1.591        6.928      4.260    0.059

From these results we can observe several things:

* in active data the time at which trajectory achieved maximum deviation
towards the wrong response seems higher in incogruent trials; however,
given that maximum deviation itself is not significant, it is unclear how to
interpret this result
* in active data also maximum y position seems significantly higher in 
incongruent trials perhaps indicating that the decision to click was easier
to reach in the congruent trials; however, given the absence of other significant
results this outcome on its own is a rather weak evidence for an effect of the trial
type
* the remaining significant differences come from passive data with small
effect sizes for maximum deviation and flips in y position, and a very small effect 
for initiation time (probably an artifact)

We can also plot the measures of interest and measures that seem to show
significant difference.


```r
N <- length(unique(mtdata$data$personid))
df <- N-1

aggregated.active <- filter(aggregated.measures, role.type=='active')

aggregated.active %>%
    group_by(personid) %>%
    select(-c(role.type, trial.type)) %>%
    summarize_all(.funs="mean") -> active.summary
grandMeans.active <- colMeans(active.summary)

adjustments.active <- active.summary
for (row in 1:dim(adjustments.active)[1]) {
    adjustments.active[row, 2:dim(adjustments.active)[2]] <- 
        grandMeans.active[2:length(grandMeans.active)] - 
        adjustments.active[row, 2:dim(adjustments.active)[2]]
}
adjustments.active <- rep(adjustments.active, 2)
adjustments.active <- arrange(adjustments.active, personid)

adjusted.active <- aggregated.active
for (row in 1:dim(adjustments.active)[1]) {
    adjusted.active[row, 4:dim(aggregated.active)[2]] <- 
        adjusted.active[row, 4:dim(aggregated.active)[2]] + 
        adjustments.active[row, 2:dim(adjustments.active)[2]]
}

adjusted.active %>%
    select(-role.type) %>%
    group_by(trial.type) %>%
    select(-personid) %>%
    summarize_all(.funs = c("mean","sd")) ->
    active.adjusted.summary

# same for passive data
aggregated.passive <- filter(aggregated.measures, role.type=='passive')
aggregated.passive %>%
    group_by(personid) %>%
    select(-c(role.type, trial.type)) %>%
    summarize_all(.funs="mean") -> passive.summary
grandMeans.passive <- colMeans(passive.summary)

adjustments.passive <- passive.summary
for (row in 1:dim(adjustments.passive)[1]) {
    adjustments.passive[row, 2:dim(adjustments.passive)[2]] <- 
        grandMeans.active[2:length(grandMeans.passive)] - 
        adjustments.passive[row, 2:dim(adjustments.passive)[2]]
}
adjustments.passive <- rep(adjustments.passive, 2)
adjustments.passive <- arrange(adjustments.passive, personid)

adjusted.passive <- aggregated.passive
for (row in 1:dim(adjustments.passive)[1]) {
    adjusted.passive[row, 4:dim(aggregated.passive)[2]] <- 
        adjusted.passive[row, 4:dim(aggregated.passive)[2]] + 
        adjustments.passive[row, 2:dim(adjustments.passive)[2]]
}

adjusted.passive %>%
    select(-role.type) %>%
    group_by(trial.type) %>%
    select(-personid) %>%
    summarize_all(.funs = c("mean","sd")) ->
    passive.adjusted.summary


compute.errorbars <- function(measure_mean, measure_sd) {
    se <- sqrt(measure_sd^2/N)
    limits.raw <- c(measure_mean[1] + c(-1,1) * qt(.975, df) * se[1],
                    measure_mean[2] + c(-1,1) * qt(.975, df) * se[2])
    return(limits.raw)
}


bars.raw <- compute.errorbars(active.adjusted.summary$AUC_mean, 
                              active.adjusted.summary$AUC_sd)
bars <- aes(ymax = bars.raw[c(2, 4)], ymin = bars.raw[c(1, 3)])
ggplot(data=active.adjusted.summary, aes(x=trial.type, AUC_mean)) +
    geom_bar(stat = "identity") +
    geom_errorbar(bars, width=0.25) +
    ggtitle("Area under curve in active data")
```

![](SocialSimon_files/figure-html/measure_plots-1.png)<!-- -->

```r
bars.raw <- compute.errorbars(active.adjusted.summary$RT_mean, 
                              active.adjusted.summary$RT_sd)
bars <- aes(ymax = bars.raw[c(2, 4)], ymin = bars.raw[c(1, 3)])
ggplot(data=active.adjusted.summary, aes(x=trial.type, RT_mean)) + 
    geom_bar(stat = "identity") +
    geom_errorbar(bars, width=0.25) +
    ggtitle("Reaction times in active data")
```

![](SocialSimon_files/figure-html/measure_plots-2.png)<!-- -->

```r
bars.raw <- compute.errorbars(passive.adjusted.summary$ypos_flips_mean, 
                              passive.adjusted.summary$ypos_flips_sd)
bars <- aes(ymax = bars.raw[c(2, 4)], ymin = bars.raw[c(1, 3)])
ggplot(data=passive.adjusted.summary, aes(x=trial.type, ypos_flips_mean)) + 
    geom_bar(stat = "identity") +
    geom_errorbar(bars, width=0.25) +
    ggtitle("Y flips in passive data")
```

![](SocialSimon_files/figure-html/measure_plots-3.png)<!-- -->


# Distributional analyses

Sometimes, averaging trajectories produces artifical results. For example, a smooth
average trajectory for a given participant could be a result of a large number
of straight trajectories that go directly to the target and discrete error type
of trajectories where participant first moves directly to the wrong side and then
abruptly changes direction (as one of the members of couple 7).

The main method to eliminate this possibility relies on bimodality analysis.
Another method includes mapping data to "trajectory prototypes".

## Bimodality analysis

This analysis checks if any of the spatial measures are bimodally distributed.


```r
# standardize measures per participant
activedata <- mt_standardize(activedata, use_variables = c("MAD", "AUC", "AD"), 
                             within = "personid", prefix = "z_")

# merge trial level data (needed for distribution qplot with facets)
mt.merged <- merge(activedata$data, activedata$measures, by="mt_id")

# plot distributions
qplot(x=z_MAD, data=mt.merged, bins=50) + facet_grid(trial.type ~ .)
```

![](SocialSimon_files/figure-html/bimodality-1.png)<!-- -->

```r
qplot(x=z_AUC, data=mt.merged, bins=50) + facet_grid(trial.type ~ .)
```

![](SocialSimon_files/figure-html/bimodality-2.png)<!-- -->

```r
qplot(x=z_AD, data=mt.merged, bins=50) + facet_grid(trial.type ~ .)
```

![](SocialSimon_files/figure-html/bimodality-3.png)<!-- -->

```r
# calculate bimodality coefficient
bm.check <- mt_check_bimodality(activedata, 
                                use_variables = c("z_MAD", "z_AUC", "z_AD"), 
                                grouping_variables = "trial.type", methods = "BC")
bm.check
```

```
## $BC
##    trial.type     z_MAD     z_AUC      z_AD
## 1   congruent 0.2540758 0.2168704 0.1862008
## 2 incongruent 0.2429527 0.2431667 0.2104477
```

A distribution is considered bimodal if BC > 0.555. In our case both distribution
plots and bimodality coefficient give no reason to suspect bimodal data that could
blur our results.


## Trajectory prototypes

Mousetrap package provides a possibility for mapping collected trajectories
to a number of trajectory prototypes frequently encountered in mouse tracking
experiments.


```r
activedata <- mt_spatialize(activedata, n_points=50)

proto <- mt_remap_symmetric(mt_prototypes, remap_xpos = "right")
mt_plot(proto, facet_col="mt_id") +
    facet_grid(.~factor(mt_id,levels=unique(mt_id))) +
    ggtitle("Trajectory prototypes")
```

![](SocialSimon_files/figure-html/prototypes-1.png)<!-- -->

```r
activedata <- mt_map(activedata, use="sp_trajectories", prototypes = proto)
mt_plot(activedata, use="sp_trajectories", 
        use2="prototyping",
        facet_col="prototype_label")
```

![](SocialSimon_files/figure-html/prototypes-2.png)<!-- -->
Even in the presence of some variability of trajectories, it seems that their
vast majority is straight.

# Trajectory dissection with PCA

Mouse trajectories are typically thought to result from a number of cognitive 
processes that depend on different factors within a trial. In order to examine
the contribution of these processes to the overall trajectory, two main methods
have been proposed: Principal Components Analysis (PCA) and visualization of beta
weights.

## PCA

To conduct PCA we use time-normalized trajectories averaged within each participant.
We will focus on components underlying x coordinates.


```r
traj.avg <- mt_aggregate_per_subject(mtdata, 
                                        use="tn_trajectories",
                                        use2_variables=c("role.type","trial.type"),
                                        subject_id="personid")
ta <- select(traj.avg, personid, role.type, trial.type, mt_seq, xpos)
ta$grp <- paste(ta$personid, ta$role.type, ta$trial.type)

ggplot(ta, aes(x=mt_seq, y=xpos, color=role.type, linetype=trial.type, group=grp)) + 
    geom_path() + ggtitle("Average trajectories of all participants")
```

![](SocialSimon_files/figure-html/pca-1.png)<!-- -->

```r
mt_plot_aggregate(mtdata, use = "tn_trajectories", 
                  x = "steps", y = "xpos", 
                  color = "trial.type", linetype="role.type",
                  subject_id = "personid") +
    theme(legend.position=c(.75,.25), legend.box='horizontal') + 
    ggtitle("Average trajectories combined by role and trial type")
```

![](SocialSimon_files/figure-html/pca-2.png)<!-- -->

```r
# transpose into wide format
traj.wide <- dcast(ta, formula=personid + role.type + trial.type ~ mt_seq, 
                   value.var="xpos")
tw <- traj.wide[,-4]

# separate by trial and role types into different data frames
tw %>% filter(., role.type=="active", trial.type=="congruent") %>%
    select(-c(personid, role.type, trial.type)) -> tw.act.congr
tw %>% filter(., role.type=="passive", trial.type=="congruent") %>%
    select(-c(personid, role.type, trial.type)) -> tw.pass.congr
tw %>% filter(., role.type=="active", trial.type=="incongruent") %>%
    select(-c(personid, role.type, trial.type)) -> tw.act.incongr
tw %>% filter(., role.type=="passive", trial.type=="incongruent") %>%
    select(-c(personid, role.type, trial.type)) -> tw.pass.incongr

# run PCA on coordinates   
pca.act.congr <- prcomp(tw.act.congr, center = TRUE, scale = TRUE)
pca.pass.congr <- prcomp(tw.pass.congr, center = TRUE, scale = TRUE)
pca.act.incongr <- prcomp(tw.act.incongr, center = TRUE, scale = TRUE)
pca.pass.incongr <- prcomp(tw.pass.incongr, center = TRUE, scale = TRUE)

var.explained <- function(sdev) {
    pr_var <- sdev^2
    prop_varex <- pr_var[1:3]/sum(pr_var) * 100
    return(prop_varex)
}

explained <- matrix(nrow=4, ncol=3)
explained[1,] <- var.explained(pca.act.congr$sdev)
explained[2,] <- var.explained(pca.pass.congr$sdev)
explained[3,] <- var.explained(pca.act.incongr$sdev)
explained[4,] <- var.explained(pca.pass.incongr$sdev)
explained <- as.data.frame(explained)
explained$grp <- c('act.congr', 'pass.congr', 'act.incongr', 'pass.incongr')

# select top 3 components from each type of trial
components <- data.frame("comp1.act.congr" = pca.act.congr$rotation[,1],
                         "comp2.act.congr" = pca.act.congr$rotation[,2],
                         "comp3.act.congr" = pca.act.congr$rotation[,3],
                         "comp1.pass.congr" = pca.pass.congr$rotation[,1],
                         "comp2.pass.congr" = pca.pass.congr$rotation[,2],
                         "comp3.pass.congr" = pca.pass.congr$rotation[,3],
                         "comp1.act.incongr" = pca.act.incongr$rotation[,1],
                         "comp2.act.incongr" = pca.act.incongr$rotation[,2],
                         "comp3.act.incongr" = pca.act.incongr$rotation[,3],
                         "comp1.pass.incongr" = pca.pass.incongr$rotation[,1],
                         "comp2.pass.incongr" = pca.pass.incongr$rotation[,2],
                         "comp3.pass.incongr" = pca.pass.incongr$rotation[,3])

components.long <- melt(components)
components.long$steps <- rep(c(1:100), 12)
components.long <- separate(components.long, variable, 
                            into = c("component", "role", "condition"))
ggplot(components.long, aes(x=steps, y=value, color=component)) + 
    geom_path() + facet_grid(role~condition) + 
    ggtitle("Top 3 principal components by trial and role type")
```

![](SocialSimon_files/figure-html/pca-3.png)<!-- -->

Given the resulting plots, we could surmise that the first component, which explains
about 77% of variance in active trials and 
83% in passive trials, reflects the constant tendency to move in a certain direction. The second component, which explains about 
14% in all types of trials decreases in congruent 
active trials while it increases in all other types. Finally, the third
component that explains 8% in active and 
2% in passive trials has also different shape between these
two types of trials.

## Beta weights

Following Scherbaum et al., we also examine trajectory angles. (TODO)
They have applied the following steps to visualize different trajectory components:

1. Angles were standardized for each participant to be between -1 and 1.
2. Created two bins of trials by a split at the median RT for each subject 
(bin 1, fast trials: M(RT) = 501 ms; bin 2, slow trials:M(RT) = 652 ms).
3. Coded four predictors for all trials: directionN(left/right), 
locationN (left/right), responseN-1(left/right), congruencyN-1. Predictors were 
coded with values -1 and 1 for easier interpretation.
4. Performed 100 multiple regressions for 4 predictors on 100 time steps separately
for each participant, which yielded 4 time-varying beta weights for each participant.
5. Computed grand average of these beta weights across participants.
6. Strength of each peak tested with one-sample t-tests of the peak beta weight
against 0.



```r
# plot trials split by 4 predictors
# plot beta weights against time
```


# Dynamical analyses

Finally, given that we are working with continuous data, we can try to apply
analyses that come from dynamical systems approach to cognition to see if
such measures can detect regularities that are missed by traditional statistical
methods.

## Sample entropy

In some cases, it may be helpful to measure the complexity of mouse trajectories. 
For example, if both response alternatives simultaneously attract participants'
mouse movement (relative to only one), this additional stress might manifest 
as less smooth, more complex, and fluctuating trajectories. 
Some mouse-tracking studies have used "x-flips", others opted for sample entropy,
a measure of predictability of trajectory given a number of surrounding coordinates.
We have seen above that a measure of x-flips did not deliver statistically significant
results. Here we calculate sample entropy, using default settings.


```r
# first for active data
# calculate sample entropy for lag 3
activedata <- mt_sample_entropy(activedata, use="tn_trajectories", m=3, 
                                dimension="xpos")
# aggregate per participant
agg.entropy.active <- mt_aggregate_per_subject(activedata, subject_id = "personid",
                                        use_variables = "sample_entropy", 
                                        use2_variables = "trial.type")
# test the difference
mt_aggregate(activedata, subject_id = "personid",
             use_variables = "sample_entropy", use2_variables = "trial.type")
```

```
##    trial.type sample_entropy
## 1   congruent     0.09562019
## 2 incongruent     0.09761750
```

```r
t.test(sample_entropy~trial.type, data=agg.entropy.active, paired=TRUE)
```

```
## 
## 	Paired t-test
## 
## data:  sample_entropy by trial.type
## t = -1.6069, df = 21, p-value = 0.123
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.0045821345  0.0005875044
## sample estimates:
## mean of the differences 
##            -0.001997315
```

```r
# for passive data
# calculate sample entropy for lag 3
passivedata <- mt_sample_entropy(passivedata, use="tn_trajectories", m=3, 
                                dimension="xpos")
# aggregate per participant
agg.entropy.passive <- mt_aggregate_per_subject(passivedata, subject_id = "personid",
                                        use_variables = "sample_entropy", 
                                        use2_variables = "trial.type")
# test the difference
mt_aggregate(passivedata, subject_id = "personid",
             use_variables = "sample_entropy", use2_variables = "trial.type")
```

```
##    trial.type sample_entropy
## 1   congruent     0.08395626
## 2 incongruent     0.08638216
```

```r
t.test(sample_entropy~trial.type, data=agg.entropy.passive, paired=TRUE)
```

```
## 
## 	Paired t-test
## 
## data:  sample_entropy by trial.type
## t = -1.6863, df = 21, p-value = 0.1065
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.0054176741  0.0005658721
## sample estimates:
## mean of the differences 
##            -0.002425901
```

```r
# TODO: entropy along y axis, angles?
```
It would seem that there is no difference in complexity of x coordinates
between different trial types in either active or passive data.


## Coupling

Another question one might ask in a joint action scenario is whether participants'
responses are more correlated within an experimental pair than across different
pairs and whether this depends on the availability of visual information about
the co-actor's movements. 

The measure that is most frequently used in the DST community is cross-recurrence
quantification analysis (CRQA) that is an index of the coupling between two time
series. Crudely put, it relies on reconstructing phase spaces of the systems from 
a given data and checking whether the states that the systems visit are close to 
each other.

In order to carry out such analysis, 3 hyperparameters need to be set or determined:

* radius: cutoff boundary that will determine if two points are recurrent or not
* delay: how many points to consider when looking for recurrence
* embedding dimension: lag unit

This can be done either based on the literature or sampled from and optimized for
a number of trials in the data and then applied to the remaining trials. Given
that we are aware of no studies that have applied CRQA to mouse-tracking data,
we opt for the latter option. We estimate the parameters based on a random sample
of trials and compute CRQA measures for all trials for all pairs.


```r
par = list(lgM =  20,
           # percentage of false neighbours excluded when adding one more 
           # dimension, relative to the first dimension
           fnnpercent = 10,
           # the span of radius examined
           # relative to SD of the distance matrix
           radiusspan = 80,       
           # the num of samples of equally spaced radius values to examine
           radiussample = 20, 
           normalize = 0, rescale = 1, mindiagline = 2, minvertline = 2,
           tw = 0, whiteline = FALSE, recpt = FALSE, typeami = "mindip")

# get random trials from a couple
get.sample.trials <- function(couple.data) {
    unique.trials <- unique(couple.data$trial)
    sample.trials <- sample(unique.trials, 20)
    return(sample.trials)
    }

exp3 %>% select(pair, personid, person, trial, x, y) %>%
    group_by(pair) %>% filter(., trial %in% get.sample.trials(.)) -> 
    sampled.data

sampled.params <- data.frame()
couples <- as.integer(unique(sampled.data$pair))

# calculate optimized parameters for sampled trials
for (num in couples) {
    couple.data <- sampled.data[sampled.data$pair==num,]
    trials <- as.integer(unique(couple.data$trial))
    for (num2 in trials) {
        trial.data <- couple.data[couple.data$trial==num2,]
        p1 <- trial.data[trial.data$person=='p1',]
        p2 <- trial.data[trial.data$person=='p2',]
        params <- try(as.data.frame(optimizeParam(p1$x, p2$x, par)), 
                      silent=TRUE)
        if ('try-error' %in% class(params)) {
            next
        }
        else if (dim(params)[1] == 0) {
            next
        }
        else {
            params$pair <- num
            params$trial <- num2
            sampled.params <- rbind(sampled.params, params)
        }
    }
}

# take the radius, delay and embedding dimensions to be the median of
# all obtained values
r = mean(sampled.params$radius)
d = median(sampled.params$delay)
e = median(sampled.params$emddim)

radius = r; delay = d; embed =  e; 
rescale =  1; normalize = 0; 
minvertline = 2; mindiagline = 2; 
whiteline = FALSE; recpt = FALSE; tw = 0


# get RQA measures from all trials
get.rqa <- function(trial.data) {
    persons <- unique(trial.data$personid)
    p1 <- trial.data[trial.data$personid==persons[1],]
    p2 <- trial.data[trial.data$personid==persons[2],]
    res <- crqa(p1$x, p2$x, delay, embed, rescale, radius,
                normalize, minvertline, mindiagline, tw,  
                whiteline, recpt)
    metrics <- data.frame(percentRecurrence = res[1], 
                 percentDeterminism = res[2],
                 longestLine = res[4],
                 entropy = res[7])
    metrics$pair <- trial.data$pair[1]
    metrics$trial <- trial.data$trial[1]
    return(metrics)
}

exp3 %>% group_by(pair, trial) %>% do(get.rqa(.)) ->
    all.metrics
```

Having computed the CRQA metrics, one can visualize the different trials, as
well as perform statistical analyses on the measures obtained.

RQA plots visualize how the states of the system evolve over time, while the
measures we can look at are as follows:

| Measure | Definition |
|:-----------------------|:-----------------------------------------------|
|recurrence rate|               how often the system visits the same state|
|determinism|                   how often the same sequences repeat|
|meanline and maxline|          how long are repeating sequences|
|entropy|                       how many repeating patterns are there|
|laminarity and trapping time|  how long the time series remains in the same state|
|trend|                         whether it's stationary|

The most widely used of these are recurrence rate (RR) and determinism.

Just to get a feel for these measures we can first examine the example trials in 
which coupling (as indexed by RR and determinism) has been estimated to be low 
and high.


```r
highestRec <- all.metrics[which(all.metrics$RR == max(all.metrics$RR)),]

hTr <- filter(exp3, pair==highestRec$pair, trial==highestRec$trial)
resH <- crqa(hTr[hTr$person=='p1',]$x, hTr[hTr$person=='p2',]$x, 
             delay, embed, rescale, radius,
             normalize, minvertline, mindiagline, tw,  whiteline, recpt)
# plot
RP <- resH$RP
RP <- matrix(as.numeric(RP), nrow = ncol(RP))
cols <- c("white","blue4")
# image(RP, xlab = "", ylab = "", col = cols)

ggplot(hTr, aes(x=x, y=y, color=person)) + geom_path() +
    xlim(0, screen.width) + ylim(0, screen.height) +
    ggtitle("Raw trajectories in a trial with highest recurrence rate")
```

![](SocialSimon_files/figure-html/rqa_plots-1.png)<!-- -->

```r
all.metrics.positive <- all.metrics[which(all.metrics$RR > 0),]
lowestRec <- all.metrics.positive[which(all.metrics.positive$RR == 
                                   min(all.metrics.positive$RR)),]

lTr <- filter(exp3, pair==lowestRec$pair, trial==lowestRec$trial)
resL <- crqa(lTr[lTr$person=='p1',]$x, lTr[lTr$person=='p2',]$x, 
             delay, embed, rescale, radius,
             normalize, minvertline, mindiagline, tw,  whiteline, recpt)

# plot
RP <- resL$RP
RP <- matrix(as.numeric(RP), nrow = ncol(RP))
cols <- c("white","blue4")
# image(RP, xlab = "", ylab = "", col = cols)

ggplot(lTr, aes(x=x, y=y, color=person)) + geom_path() +
    xlim(0, screen.width) + ylim(0, screen.height) +
    ggtitle("Raw trajectories in a trial with lowest recurrence rate")
```

![](SocialSimon_files/figure-html/rqa_plots-2.png)<!-- -->

```r
highestDET <- all.metrics[which(all.metrics$DET == max(all.metrics$DET,
                                                       na.rm=TRUE)),]

hTrDET <- filter(exp3, pair==highestDET$pair[1], trial==highestDET$trial[1])
resHDET <- crqa(hTrDET[hTrDET$person=='p1',]$x, hTrDET[hTrDET$person=='p2',]$x, 
             delay, embed, rescale, radius,
             normalize, minvertline, mindiagline, tw,  whiteline, recpt)
# plot
RP <- resHDET$RP
RP <- matrix(as.numeric(RP), nrow = ncol(RP))
cols <- c("white","blue4")
# image(RP, xlab = "", ylab = "", col = cols)

ggplot(hTrDET, aes(x=x, y=y, color=person)) + geom_path() +
    xlim(0, screen.width) + ylim(0, screen.height) +
    ggtitle("Raw trajectories in a trial with highest determinism")
```

![](SocialSimon_files/figure-html/rqa_plots-3.png)<!-- -->

```r
lowestDET <- all.metrics.positive[which(all.metrics.positive$DET == 
                                            min(all.metrics.positive$DET, 
                                                na.rm=TRUE)),]

lTrDET <- filter(exp3, pair==lowestDET$pair, trial==lowestDET$trial)
resL <- crqa(lTrDET[lTrDET$person=='p1',]$x, lTrDET[lTrDET$person=='p2',]$x, 
             delay, embed, rescale, radius,
             normalize, minvertline, mindiagline, tw,  whiteline, recpt)
# plot
RP <- resL$RP
RP <- matrix(as.numeric(RP), nrow = ncol(RP))
cols <- c("white","blue4")
# image(RP, xlab = "", ylab = "", col = cols)

ggplot(lTrDET, aes(x=x, y=y, color=person)) + geom_path() +
    xlim(0, screen.width) + ylim(0, screen.height) +
    ggtitle("Raw trajectories in a trial with lowest determinism")
```

![](SocialSimon_files/figure-html/rqa_plots-4.png)<!-- -->

Subsequently we can see if coupling is higher within couples than across and whether
it is affected by the presence of visual information (TODO: compare to coupling
in condition 4).

Judging whether coupling within couples is indeed present typically implies a 
comparison between coupling calculated for real couples and coupling calculated
for so-called surrogate (virtual, fake) couples. That is, we take data from the
same people but form new pairs from them and compute the same CRQA measures. Of
course, given that people perform the same task, some amount of coupling is expected
to hold just by virtue of the task requirements. However, if there is significant
coupling between people who perform the task together, it should be higher between
people who actually performed it together than between people who merely performed
the same task.


```r
all.metrics <- ungroup(all.metrics)
all.metrics %>% 
    group_by(pair) %>% 
    summarize_all(funs(mean(., na.rm=TRUE))) %>%
    select(-trial) ->
    all.metrics.summary
all.metrics.summary$cond <- "real"

exp3 %>%
  group_by(personid) %>%
  mutate(count = n()) ->
    exp3

exp3$pair <- as.integer(as.character(exp3$pair))
pairids <- rep(unique(exp3$pair), each=2)
counts <- rep(unique(exp3$count), each=2)

i <- 0
num.tests <- 10
while (i < num.tests) {
    df <- data.frame(count=counts, pairs=sample(pairids))
    newpairs <- df[rep(seq_len(dim(df)[1]), df$count), 2]
    exp3$fake <- newpairs
    
    exp3 %>% group_by(fake, trial) %>% 
        do(get.rqa(.)) -> 
        all.metrics
    all.metrics <- ungroup(all.metrics)
    
    all.metrics %>% 
        group_by(pair) %>% 
        summarize_all(funs(mean(., na.rm=TRUE))) %>%
        select(-c(trial, fake)) ->
        new.metrics.summary
    new.metrics.summary$cond <- paste0("fake", i)
    all.metrics.summary <- rbind(all.metrics.summary, new.metrics.summary)
    
    i <- i+1
}

all.metrics.summary$cond <- as.factor(all.metrics.summary$cond)

outcomes <- matrix(nrow=num.tests, ncol = 2)
for (i in 0:(num.tests-1)) {
    surrogate <- paste0('fake',i)
    dat <- filter(all.metrics.summary, cond==surrogate | cond=='real')
    res <- t.test(data=dat, RR~cond, paired=FALSE)
    new.row <- c(surrogate, round(res$p.value, 5))
    outcomes[i+1,] <- new.row
}
outcomes <- as.data.frame(outcomes)
colnames(outcomes) <- c("comparison.group", "p.value")
outcomes
```

```
##    comparison.group p.value
## 1             fake0 0.90845
## 2             fake1 0.33642
## 3             fake2 0.64938
## 4             fake3 0.90299
## 5             fake4 0.72614
## 6             fake5 0.31054
## 7             fake6 0.47386
## 8             fake7 0.53362
## 9             fake8   0.943
## 10            fake9 0.95376
```

```r
ggplot(all.metrics.summary, aes(x=cond, y=RR)) + 
    geom_bar(stat = "summary", fun.y = "mean") +
    ggtitle("Recurrence rate in real and fake couples")
```

![](SocialSimon_files/figure-html/rqa_test-1.png)<!-- -->

From the plot as well as t-test it seems that real couples do not show higher
level of coupling than fake pairs. This, together with a lack of social Simon
effect in the standard trajectory measures, as well as trajectories showing spatial
division of labor, suggests that people performed the task individually, rather 
than approaching it as a social, joint activity.

## Fractal analysis

Fractal measures detrmine how variability scales with sample size. Low level
of self-similarity typically indicates a random process, medium level - 
self-organization, while high level an influence of external constraints on
the unfolding of cognitive dynamics.
(TODO)

