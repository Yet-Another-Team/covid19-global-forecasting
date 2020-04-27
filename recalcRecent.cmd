timeout 12600
cd C:\ML\COVID-19\
git pull
cd C:\ML\
copy /Y C:\ML\COVID-19\csse_covid_19_data\csse_covid_19_daily_reports C:\ML\covid19-global-forecasting\jhu_data\csse_covid_19_daily_reports
copy /Y C:\ML\COVID-19\csse_covid_19_data\csse_covid_19_time_series C:\ML\covid19-global-forecasting\jhu_data\csse_covid_19_time_series
cd C:\ML\covid19-global-forecasting
alpheus compute 2_SIR_spline_beta\params.csv