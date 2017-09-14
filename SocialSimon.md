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
exp1 <- load.experiment(1)
exp2 <- load.experiment(2)
exp3 <- load.experiment(3)
exp4 <- load.experiment(4)
```

Then we add the coding for independent variables on 

* trial type, which depends on the congruency between the color of the cue and 
its location
* whose turn it was to respond in a given trial, which depends on the color of the cue
* previous trial type and whose turn it was to respond




```r
exp1 <- add.current.vars(exp1)
exp2 <- add.current.vars(exp2)
exp3 <- add.current.vars(exp3)
exp4 <- add.current.vars(exp4)
```




```r
exp1 <- add.previous.vars(exp1)
exp2 <- add.previous.vars(exp2)
exp3 <- add.previous.vars(exp3)
exp4 <- add.previous.vars(exp4)
```

The y-coordinates are immediately flipped vertically because the package that was
used for collecting the data (Matlab Psychtoolbox) encodes the screen's top left 
as coordinates [0, 0] and therefore y-coordinates grow towards the bottom of the screen
while for ease of analysis we would like them to grow towards the screen's top.

We plot an example trial before any further pre-processing. The complete trajectory 
that is plotted contains all coordinates recorded since the start of the
trial and so also before the participants clicked on the start button located
at the bottom of the screen.

![](SocialSimon_files/figure-html/07_flip_y-1.png)<!-- -->

Now we filter out only successful trials, in which participants did not miss any
deadlines and a correct response was given.




```r
exp1.complete <- remove.failed(exp1)
exp2.complete <- remove.failed(exp2)
exp3.complete <- remove.failed(exp3)
exp4.complete <- remove.failed(exp4)
```

Complete trajectories will be needed for the dynamical part of our analysis at 
the end of this script. For remaining analyses we need to extract
the portions of the trajectories after participants have clicked on the start
button.




```r
exp1.subset <- get.subset(exp1.complete)
exp2.subset <- get.subset(exp2.complete)
exp3.subset <- get.subset(exp3.complete)
exp4.subset <- get.subset(exp4.complete)
```

![](SocialSimon_files/figure-html/12_plot_extracted-1.png)<!-- -->

For convenience we rescale the coordinates into a standard MouseTracker 
coordinate space, where x is in range [-1, 1] and y in range [0, 1.5].




```r
exp1.scaled <- get.scaled(exp1.subset)
exp2.scaled <- get.scaled(exp2.subset)
exp3.scaled <- get.scaled(exp3.subset)
exp4.scaled <- get.scaled(exp4.subset)
```

![](SocialSimon_files/figure-html/15_plot_rescaled-1.png)<!-- -->

Now we align all the trajectories to the common [0, 0] origin and timestamps
to start at 0.




```r
exp1.aligned <- get.aligned(exp1.scaled)
exp2.aligned <- get.aligned(exp2.scaled)
exp3.aligned <- get.aligned(exp3.scaled)
exp4.aligned <- get.aligned(exp4.scaled)
```

![](SocialSimon_files/figure-html/18_plot_aligned-1.png)<!-- -->

At this point we can already visualize all trajectories of all individual participants
and pairs.

![](SocialSimon_files/figure-html/19_all_trajectories-1.png)<!-- -->![](SocialSimon_files/figure-html/19_all_trajectories-2.png)<!-- -->![](SocialSimon_files/figure-html/19_all_trajectories-3.png)<!-- -->![](SocialSimon_files/figure-html/19_all_trajectories-4.png)<!-- -->

We see that most pairs seem to divide the screen space between each other
by moving mostly directly towards their assigned response box and avoiding the center.
However, there are exceptions. 

In condition 3 (with visual feedback) one particular pair (number 7) has mostly 
upward moving trajectories. We can further explore here whether the joint upward 
motion is induced by one of the participants or happens immediately on both sides.

![](SocialSimon_files/figure-html/20_aberrant_pair-1.png)<!-- -->

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

In condition 4 (without visual feedback) it is rather than certain individuals 
adopt the "move upward" strategy, independently of their partner (which is to be
expected since they do not see the partner's movements).

As the next pre-processing step we will flip all trajectories to one side.
This ensures that every trajectory starts at the bottom of the coordinate system 
and ends in the top right corner. It is done to obtain comparable trajectory 
measures.




```r
exp1.clean <- get.clean(exp1.aligned)
exp2.clean <- get.clean(exp2.aligned)
exp3.clean <- get.clean(exp3.aligned)
exp4.clean <- get.clean(exp4.aligned)
```

As a last step, a look at sampling rate distribution to see whether it reveals 
any outliers that could indicate missing data or wrong recording.



From the data we can extract all intervals between adjacent sampling points and
calculate their mean and standard devation. We then establish a cutoff point of
3 SD from the mean beyond which the sampling interval is considered to be an outlier.
We calculate the number of such outliers for each condition and their mean sampling
rate. The resulting numbers are presented in the following table.


              rate.means   rate.sds   rate.cutoff   num.outliers   mean.outliers
-----------  -----------  ---------  ------------  -------------  --------------
condition1         10.87       0.01         10.88            596           10.90
condition2         10.87       0.01         10.91              8           13.60
condition3         10.87       1.30         14.78             41          207.71
condition4         10.87       0.53         12.47             39           42.92

It turns out that the means sampling rate for all conditions is 
within expected parameters given the set sampling rate of 92 Hz. The number of
outliers and their means varies per condition with much larger deviations in social
conditions. Given the cross-computer data stream in these conditions, some amount
of data loss was to be expected. In order to facilitate binned data analysis, 
we remove trials that contain these large deviations (separately for each condition).

Next we look at y coordinates that should not be lower than some margin around 0
(after flipping and alignment). A y-coordinate that is more negative indicates a
faulty recording of the start button press.





With the sampling timing issues fixed, we can examine reaction time outliers and
add another variable to the data, which indicates whether the trial was fast or slow.




```r
exp1.clean <- time.outliers(exp1.clean)
exp2.clean <- time.outliers(exp2.clean)
exp3.clean <- time.outliers(exp3.clean)
exp4.clean <- time.outliers(exp4.clean)
```

![](SocialSimon_files/figure-html/29_plot_rts-1.png)<!-- -->

Note that in condition 2 the histogram looks markedly different for active and
passive trials. This is because when participant was supposed to refrain from 
responding, they had to simply wait for the trial to end.

After we produced clean data for all conditions, we combine the two social ones 
into a single data frame.


```r
exp.clean <- rbind(exp3.clean, exp4.clean)
```

For the analysis we will use a mousetrap package. Accordingly, the next step is 
to transform the data into a mousetrap object.


```r
mtdata1 <- mt_import_long(exp1.clean, 
                        xpos_label = 'x.flipped', 
                        ypos_label = 'y.aligned', 
                        timestamps_label = 't.aligned', 
                        mt_id_label = c('personid','trial'), 
                        mt_seq_label = 'sampling.order', 
                        reset_timestamps = FALSE)
mtdata2 <- mt_import_long(exp2.clean, 
                        xpos_label = 'x.flipped', 
                        ypos_label = 'y.aligned', 
                        timestamps_label = 't.aligned', 
                        mt_id_label = c('personid','trial'), 
                        mt_seq_label = 'sampling.order', 
                        reset_timestamps = FALSE)

mtdata <- mt_import_long(exp.clean, 
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
mtdata1 <- mt_time_normalize(mtdata1)
mtdata2 <- mt_time_normalize(mtdata2)
mtdata <- mt_time_normalize(mtdata)
```

We need to check that all participants have a balanced number of observations 
for different variables of interest.




```r
counts1 <- get.counts(exp1.complete)
counts2 <- get.counts(exp2.complete)
counts3 <- get.counts(exp3.complete)
counts4 <- get.counts(exp4.complete)
```

![](SocialSimon_files/figure-html/35_plot_counts-1.png)<!-- -->![](SocialSimon_files/figure-html/35_plot_counts-2.png)<!-- -->![](SocialSimon_files/figure-html/35_plot_counts-3.png)<!-- -->![](SocialSimon_files/figure-html/35_plot_counts-4.png)<!-- -->

From the counts plot, we can see that in condition 1 and 2 one person in each
has markedly less successful trials than other participants. In condition 3, 
even though different numbers of trials remain for different pairs after removing 
unsuccessful ones, there are similar counts for different trial and role types.

By contrast, in condition 4, three of the pairs completely lose observations for
one of the incongruent set of trials, namely participants with ids 14, 23, 29 
lose incongruent active trials while their co-actors (ids 13, 24, 30) incongruent
passive ones. Further investigation of these counts reveals that the crucial step
that leads to this loss is filtering out trials in which incorrect response was
given, that is, participant clicked on the wrong response box. Given that this 
correlates with trials being incongruent (the cue appeared on the same side as
the incorrect response box), we might infer that participants 14, 23 and 29
misunderstood the task, i.e. they were responding to the location of the cue, rather
than its color. We can confirm this conclusion by plotting the trajectories of
these participants.

![](SocialSimon_files/figure-html/36_plot_wrong_people-1.png)<!-- -->

We see that indeed, 3 participants in condition 4 misunderstood the instructions
and therefore need to be removed from further analysis.


```r
mtdata1 <- mt_subset(mtdata1, personid != "c1p5")
mtdata2 <- mt_subset(mtdata2, personid != "c2p11")
mtdata2 <- mt_subset(mtdata2, !(personid == "c2p16" & role.type == "passive"))
# redo this later
mtdata2 <- mt_subset(mtdata2, !(personid == "c2p17" & role.type == "passive"))

mtdata <- mt_subset(mtdata, !(condition == 4 & (personid %in% 
                                                    c("c4p4", "c4p13", "c4p19"))))
mtdata <- mt_subset(mtdata, !(condition == 4 & role.type == 'passive' & 
                                    (personid %in% c("c4p3", "c4p14", "c4p20"))))
```




As a result of data cleaning, 13%, 3%, 
4%, 27% of trials in conditions 1 to 4 are removed.

# Dependent variables calculation

There is a number of measures that can be calculated on the basis of raw time
and normalized trajectories. 
First, we retrieve trajectory derivatives (velocity, acceleration) and angles 
from the raw time data. Next, we calculate a variety of measures on normalized
data. 


```r
# get mt measures for all experiments
mtdata1 <- mt_derivatives(mtdata1, use="trajectories")
mtdata1 <- mt_angles(mtdata1, use="trajectories")
mtdata1 <- mt_measures(mtdata1, use="trajectories", save_as = "measures")

mtdata2 <- mt_derivatives(mtdata2, use="trajectories")
mtdata2 <- mt_angles(mtdata2, use="trajectories")
mtdata2 <- mt_measures(mtdata2, use="trajectories", save_as = "measures")

mtdata <- mt_derivatives(mtdata, use="trajectories")
mtdata <- mt_angles(mtdata, use="trajectories")
mtdata <- mt_measures(mtdata, use="trajectories", save_as = "measures")
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
# we split only the social conditions by role
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

![](SocialSimon_files/figure-html/41_explo_plots_12-1.png)<!-- -->

Now for the social conditions.

![](SocialSimon_files/figure-html/42_explo_plots-1.png)<!-- -->

From these plots it would seem that neither of the independent variables
affects the shape of movement trajectory and furthermore that participants go
straight for their assigned response button instead of producing curved trajectories
typically found in mouse-tracking studies. 

However, to check whether this overall pattern holds on the individual level, 
we can also plot average trajectories for individual participants, focusing here 
on the effect of current trial type, for active and passive data.





![](SocialSimon_files/figure-html/45_avg_person-1.png)<!-- -->![](SocialSimon_files/figure-html/45_avg_person-2.png)<!-- -->![](SocialSimon_files/figure-html/45_avg_person-3.png)<!-- -->![](SocialSimon_files/figure-html/45_avg_person-4.png)<!-- -->

Here we can make several observations. First, even though the general pattern that
we see on the group level holds for majority but not all individuals. There are
3 participants in condition 3 and 2 participants in condition 4 do exhibit trajectories
that go upwards first and then to the response box (curved trajectories). For 
some of these individuals we can also note a slightly bigger curve in incogruent
active trials but no difference in passive trials. Furthermore, all of them show
shorter trajectories in passive than active trials, while most participants with 
straight trajectories seem to proceed all the way to the response box even when
it is not their turn to respond, i.e. despite the fact that they could return to 
the starting position as soon as it became clear it is not their trial.

We can confirm the latter observation by also plotting a histogram of maximum y 
coordinate reached in passive trials. The thresholds indicated on the plot are the 
locations of the start boundary, the y-coordinated that had to be crossed in order 
for the cue to appear and the lower response box boundary. As can be seen from the plot,
the majority of trajectories goes beyond that last threshold.

![](SocialSimon_files/figure-html/46_passive_hist-1.png)<!-- -->

For the following plots we will separate straight-trajectory participants from the
curved-trajectory participants in conditions 2, 3 and 4 as they seem to behave in a qualitatively different manner.



In addition to plotting time-normalized trajectories in full, we can bin them into
several intervals with the aim of further subjecting the bins to inferential testing.

![](SocialSimon_files/figure-html/48_binned_tnplots-1.png)<!-- -->![](SocialSimon_files/figure-html/48_binned_tnplots-2.png)<!-- -->![](SocialSimon_files/figure-html/48_binned_tnplots-3.png)<!-- -->![](SocialSimon_files/figure-html/48_binned_tnplots-4.png)<!-- -->![](SocialSimon_files/figure-html/48_binned_tnplots-5.png)<!-- -->

Some analyses, such as looking at movement velocity profiles, require retaining 
trajectories in raw time. In this case, we decide how many raw time bins to create 
between 0 ms and some cutoff (e.g., 1500 ms) and then create a number of raw time steps. 
Thus, each step (i.e., coordinate pair) of a trajectory reflects the location of 
the mouse during some raw time bin (e.g., 500-600 ms if bins are 100 ms wide).



Once we have created such raw time bins, we can look at velocity in different
types of trials. The literature suggests, for example, that stronger competition 
between response options should be characterized by an initial decreased velocity 
as competing choices inhibit each other, followed by an increase in velocity once 
the system converges upon a decision and the inhibition is alleviated. 
Thus, analyzing velocity data can allow for inferences about when commitments to 
a particular response are made.

In our particular case, we might also ask whether velocity is different between
active and passive trials and whether participants that exhibit qualitatively
different movement trajectories also differ in their velocity profiles.



Here we plot velocity profiles in binned raw time, together with approximate
time in which cue appeared (mean appearance time being 294 ms).

![](SocialSimon_files/figure-html/plot_velocity-1.png)<!-- -->![](SocialSimon_files/figure-html/plot_velocity-2.png)<!-- -->![](SocialSimon_files/figure-html/plot_velocity-3.png)<!-- -->![](SocialSimon_files/figure-html/plot_velocity-4.png)<!-- -->![](SocialSimon_files/figure-html/plot_velocity-5.png)<!-- -->![](SocialSimon_files/figure-html/plot_velocity-6.png)<!-- -->![](SocialSimon_files/figure-html/plot_velocity-7.png)<!-- -->![](SocialSimon_files/figure-html/plot_velocity-8.png)<!-- -->![](SocialSimon_files/figure-html/plot_velocity-9.png)<!-- -->




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


```
## `.cols` has been renamed and is deprecated, please use `.vars`
## `.cols` has been renamed and is deprecated, please use `.vars`
```

The test revealed a sequence of 0 and 4 significant t-tests in
conditions 3 and 4 respectively, on the difference between x-coordinates in congruent 
and incongruent trials. In order to determine what is the minimum number of 
significant t-tests that qualifies as a pattern, a bootstrapping procedure would 
be required. However, in previous research 8 was the minimum so this number seems 
rather low by comparison.

## Anovas on binned trajectories

Other than looking at particular coordinates, we can also run tests on binned
trajectories (both normalized and raw time), that we have also plotted above.
In this case the analysis we perform is repeated measures 3 (bins) by 2 (trial
type) ANOVA. 


```r
avg.tn.trajectory.bins <-  
    mt_aggregate_per_subject(active.straight, 
                            use="av_tn_trajectories",
                            use2_variables=c("condition", "trial.type"),
                            subject_id="personid")
avg.tn.trajectory.bins$bin <- factor(avg.tn.trajectory.bins$steps)

aov_ez(data=avg.tn.trajectory.bins,
       id="personid", dv="xpos", within = c("trial.type", "bin"),
       between="condition",
       anova_table = list(es=c("ges", "pes"), correction=c("GG")))
```

```
## Warning: Missing values for following ID(s):
## c4p13, c4p19, c4p4
## Removing those cases from the analysis.
```

```
## Warning in aov(formula(paste(dv.escaped, "~", paste(c(between.escaped,
## within.escaped), : Error() model is singular
```

```
## Anova Table (Type 3 tests)
## 
## Response: xpos
##                     Effect          df  MSE           F    ges  pes
## 1                condition       1, 32 0.02        0.04  .0008 .001
## 2               trial.type       1, 32 0.00        0.46 <.0001  .01
## 3     condition:trial.type       1, 32 0.00        2.24  .0003  .07
## 4                      bin 1.91, 61.02 0.00 1520.11 ***    .93  .98
## 5            condition:bin 1.91, 61.02 0.00        1.81    .02  .05
## 6           trial.type:bin 1.96, 62.73 0.00        0.73 <.0001  .02
## 7 condition:trial.type:bin 1.96, 62.73 0.00        0.18 <.0001 .006
##   p.value
## 1     .85
## 2     .50
## 3     .14
## 4  <.0001
## 5     .17
## 6     .48
## 7     .83
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

aov_ez(data=avg.fast.trajectory.bins, id="personid", dv="xpos", 
       within = c("trial.type","bin"),
       anova_table = list(es=c("ges","pes"), correction=c("GG")))
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




measure           type            t   df         p   lower.conf   upper.conf   estimate   effect
----------------  --------  -------  ---  --------  -----------  -----------  ---------  -------
ypos_max          active     -2.141   21   0.04416       -0.006        0.000     -0.003    0.423
MD_above          passive    -2.260   21   0.03456       -0.015       -0.001     -0.008    0.289
MD_above_time     active     -2.929   21   0.00802      -33.459       -5.675    -19.567    0.539
ypos_flips        passive    -2.186   21   0.04025       -0.090       -0.002     -0.046    0.280
initiation_time   passive     3.314   21   0.00330        1.583        6.917      4.250    0.059



measure        type            t   df         p   lower.conf   upper.conf   estimate   effect
-------------  --------  -------  ---  --------  -----------  -----------  ---------  -------
xpos_max       passive    -2.886   13   0.01275       -0.022       -0.003     -0.012    0.280
MAD_time       passive     2.357   13   0.03478        0.952       21.893     11.423    0.288
vel_max_time   passive     3.080   13   0.00878        2.763       15.744      9.253    0.457

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

bars.raw <- compute.errorbars(active.adjusted.summary$RT_mean, 
                              active.adjusted.summary$RT_sd)
bars <- aes(ymax = bars.raw[c(2, 4)], ymin = bars.raw[c(1, 3)])
ggplot(data=active.adjusted.summary, aes(x=trial.type, RT_mean)) + 
    geom_bar(stat = "identity") +
    geom_errorbar(bars, width=0.25) +
    ggtitle("Reaction times in active data")

bars.raw <- compute.errorbars(passive.adjusted.summary$ypos_flips_mean, 
                              passive.adjusted.summary$ypos_flips_sd)
bars <- aes(ymax = bars.raw[c(2, 4)], ymin = bars.raw[c(1, 3)])
ggplot(data=passive.adjusted.summary, aes(x=trial.type, ypos_flips_mean)) + 
    geom_bar(stat = "identity") +
    geom_errorbar(bars, width=0.25) +
    ggtitle("Y flips in passive data")
```


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
## 1   congruent 0.2600469 0.2264174 0.1967166
## 2 incongruent 0.2598736 0.2337273 0.2051982
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
about 69% of variance in active trials and 
75% in passive trials, reflects the constant tendency to move in a certain direction. The second component, which explains about 
20% in all types of trials decreases in congruent 
active trials while it increases in all other types. Finally, the third
component that explains 9% in active and 
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
## 1   congruent     0.09608580
## 2 incongruent     0.09717253
```

```r
t.test(sample_entropy~trial.type, data=agg.entropy.active, paired=TRUE)
```

```
## 
## 	Paired t-test
## 
## data:  sample_entropy by trial.type
## t = -1.1591, df = 38, p-value = 0.2536
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.0029847149  0.0008112553
## sample estimates:
## mean of the differences 
##             -0.00108673
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
## 1   congruent     0.08248634
## 2 incongruent     0.08433156
```

```r
t.test(sample_entropy~trial.type, data=agg.entropy.passive, paired=TRUE)
```

```
## 
## 	Paired t-test
## 
## data:  sample_entropy by trial.type
## t = -1.7216, df = 35, p-value = 0.09398
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.004021170  0.000330721
## sample estimates:
## mean of the differences 
##            -0.001845224
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
    sampled.data3

exp4 %>% select(pair, personid, person, trial, x, y) %>%
    filter(!(pair %in% c(2, 7, 10))) %>%
    group_by(pair) %>% filter(., trial %in% get.sample.trials(.)) -> 
    sampled.data4

# calculate optimized parameters for sampled trials
calculate.params <- function(sampled.data) {
    sampled.params <- data.frame()
    couples <- as.integer(unique(sampled.data$pair))
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
    return(sampled.params)
}

sampled.params3 <- calculate.params(sampled.data3)
sampled.params4 <- calculate.params(sampled.data4)

# get RQA measures from all trials
get.rqa <- function(trial.data, sampled.params) {
    # take the radius, delay and embedding dimensions to be the median of
    # all obtained values
    radius = mean(sampled.params$radius)
    delay = median(sampled.params$delay)
    embed = 2 # this most frequently comes out as 2

    rescale =  1; normalize = 0; 
    minvertline = 2; mindiagline = 2; 
    whiteline = FALSE; recpt = FALSE; tw = 0

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

exp3 %>% group_by(pair, trial) %>% do(get.rqa(., sampled.params3)) ->
    all.metrics3
exp4 %>% filter(!(pair %in% c(2, 7, 10))) %>%
    group_by(pair, trial) %>% do(get.rqa(., sampled.params4)) ->
    all.metrics4
```


```r
# test
all.metrics3 <- ungroup(all.metrics3)
all.metrics3 %>% 
    group_by(pair) %>% 
    summarize_all(funs(mean(., na.rm=TRUE))) %>%
    select(-trial) ->
    all.metrics.summary3
all.metrics.summary3$cond <- as.character(3)

all.metrics4 <- ungroup(all.metrics4)
all.metrics4 %>% 
    group_by(pair) %>% 
    summarize_all(funs(mean(., na.rm=TRUE))) %>%
    select(-trial) ->
    all.metrics.summary4
all.metrics.summary4$cond <- as.character(4)

all.metrics.summary <- rbind(all.metrics.summary3, all.metrics.summary4)

p1 <- ggplot(all.metrics.summary, aes(x=cond, y=RR)) + 
    geom_bar(stat = "summary", fun.y = "mean") +
    ggtitle("Recurrence rate in two social conditions")
p2 <- ggplot(all.metrics.summary, aes(x=cond, y=DET)) + 
    geom_bar(stat = "summary", fun.y = "mean") +
    ggtitle("Recurrence rate in two social conditions")
p3 <- ggplot(all.metrics.summary, aes(x=cond, y=maxL)) + 
    geom_bar(stat = "summary", fun.y = "mean") +
    ggtitle("Recurrence rate in two social conditions")
p4 <- ggplot(all.metrics.summary, aes(x=cond, y=rENTR)) + 
    geom_bar(stat = "summary", fun.y = "mean") +
    ggtitle("Recurrence rate in two social conditions")

multiplot(p1, p2, p3, p4, cols=2)
```

![](SocialSimon_files/figure-html/test_social_rqas-1.png)<!-- -->

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

