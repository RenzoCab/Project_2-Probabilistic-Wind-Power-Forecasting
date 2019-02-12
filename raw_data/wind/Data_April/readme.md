From: Raul Tempone <rtempone@gmail.com>
Subject: new wind data from Uruguay with proper normalization and adjustment of the initial point
Date: April 3, 2017 at 6:44:50 PM GMT+3
To: Evangelia Kalligiannaki <evangelia.kalligiannaki@kaust.edu.sa>, soumaya elkantassi <soumaya.elkantassi@gmail.com>, Soumaya Elkantassi <soumaya.kantassi@kaust.edu.sa>

Dear both,


In the attached zip file you can find the data that we were promised.

The important files are four:

Sal_gahas_h72_2016_130317_norm_asim.txt
Sal_mtlog_h72_2016_130317_norm_asim.txt
Sal_utep3_h72_2016_130317_norm_asim.txt
Sal_uter1_2016_130317_iniper_norm.txt

They correspond to the forecasts of Garrand Hassam, Meteorological, UTE / IMFIA and the actual data of wind power. For each day four forecasts are generated at 72 hours, one every six hours, specifically at: 01:00, 7:00, 13:00 and 19:00.

Our uruguayan contact has left the data cleansed in two senses. First, they are normalized by the total installed power, i.e., the power data is always expressed as a fraction [0,1] between power and installed capacity. Whenever possible, they ran the power assimilation algorithm, to "paste" the forecast time series, to the initial condition. I understand that it was not possible between September 15 and 20, 2016.

In addition to the forecast itself (POTP), companies add a band of confidence. We do not know how Garrand Hassam and Meteorological calculate it. In the case of the IMFIA, it arises from a sensitivity analysis in the initial condition (disturbing the initial condition, obtaining 60 trajectories, and taking the upper and lower envelopes). I imagine that this data will not be used in our computations, but they left it because it is interesting to compare against it. .

Please let me know if you have questions.


All the best
Raul
-- 
------------------------------
R. Tempone, PhD
Professor of Applied Mathematics;
Division of Mathematics & Computational Sciences & Engineering (CEMSE)
Building #1 (UN 1550), Office No 4109, Level 4
4700 King Abdullah University of Science and Technology
Thuwal 23955-6900, Kingdom of Saudi Arabia
http://www.stochastic_numerics.kaust.edu.sa
http://www.kaust.edu.sa
------------------------------
