## Assassination Classroom Exercise

#### Installation

This project uses the following dependencies in R.

* tidyverse
* ggplot2
* jsonlite
* sp
* rgdal (likely optional)
* dlm

#### Running

The prediction logic is located in `bomb_plan.R`, but to execute, the `main.R`
file provides utilities to easily run prediction given an array of start times
and start locations.

The steps to execute are thus:

1. Open `main.R`, and at the top, modify the `start_time`, `start_latitude` and
   `start_longitude` vectors with your data. `start_time` values should be of the form, 
   `"2020-08-18T17:50:42.000Z"` and `start_latitude` and `start_longitude`
   should be floating point values of standard latitude and longitude values.
2. After adding data as above, run `Rscript main.R` and it will print to console a
   prediction in the following form.
   ```
   Day Index: 1

   PREDICTION 1:
   Location (lat,long):
   UTM Location (lat,long):
   UTC Time (seconds since epoch):

   PREDICTION 2:
   Location (lat,long):
   UTM Location (lat,long):
   UTC Time (seconds since epoch):
   ```
3. Alternatively, you can programmatically access the prediction by importing
   the `easy_predict` function from `bomb_plan.R` and then passing it three
   vectors of the form described in step 1. It will return a list of the form:
   ```
   list(bomb1=list(utm_latitude=..., utm_longitude=...,
                   latitude=..., longitude=...,
                   time=...),
        bomb2=list(utm_latitude=..., utm_longitude=...,
                   latitude=..., longitude=...,
                   time=...))
   ```

## Performance

On my own computer, the script takes approximately 6 seconds to produce a
prediction.
