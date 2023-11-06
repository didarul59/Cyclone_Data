import xarray as xr
import numpy as np
import pandas as pd

from google.colab import drive
drive.mount('/content/drive')

import xarray as xr

# Open the NetCDF file
ds = xr.open_dataset('/content/drive/MyDrive/data/adaptor.mars.internal-1694757024.7621102-8867-12-d7dc3167-39d3-47b1-b6a4-c69e86845ad0.nc')

# Check if the "expver" dimension exists
if 'expver' in ds.dims:
    # Select a specific value of the "expver" dimension (e.g., 1)
    ds = ds.isel(expver=0)  # Replace 0 with the index of the value you want to keep

# Save the modified dataset to a new NetCDF file
ds.to_netcdf('output_file5.nc')

era_data = xr.open_mfdataset('/content/output_file5.nc').sel(latitude=slice(25, 0), longitude=slice(79,101)) #see main data
era_data

# selecting variable for lat lon and time
lon = era_data['longitude'][:]
lat = era_data['latitude'][:]
time = era_data['time'][:]

#converting lon lat in lon2d and lat2d
lon2d, lat2d = np.meshgrid(lon, lat)

u10 = era_data['u10']
v10 = era_data['v10']

wind=np.sqrt(u10**2 + v10**2)
wind

nlon = len(wind.longitude)
nlat = len(wind.latitude)
nt = len(wind.time)

wind.data

wind_data = np.reshape(wind.data,((nt, nlat*nlon)))
wind_data = wind_data.compute()
wind_data

#calculating the maximum speed and the index of the cyclone
cyclonic_wind_speed = []
index = []

for i in range(len(wind.time)):
  maximum = np.nanmax(wind_data[i,:])
  if maximum>17:
    index.append(i)
    cyclonic_wind_speed.append(maximum)

print(index)

len(cyclonic_wind_speed)

#time
datetimeindex = era_data.time.data
cyclone_date = []

for i in range(len(index)):
  cyclone_date.append(datetimeindex[index[i]])

#data_table0 is used for dummy database
#so that main database wont affect
list0 = {'cyclonic_wind_speed': cyclonic_wind_speed, 'Date': cyclone_date}
data_table0 = pd.DataFrame(list0)
data_table0.max()

#Starting and Ending date from the wind
primary_count = 0 #primary count
start = cyclone_date[0]
starting_date = []
ending_date = []


for i in range(len(cyclone_date)):
  time_difference = (cyclone_date[i] - cyclone_date[i-1]).astype('timedelta64[h]')

  if time_difference >= np.timedelta64(48, 'h'):
    primary_count = primary_count+1
    starting_date.append(start)
    ending_date.append(cyclone_date[i-1])
    start = cyclone_date[i]

print("Primary cyclone count: ", primary_count)
print(starting_date)
print(ending_date)

list_value = {'Starting Date': starting_date, 'Ending Date': ending_date}
data_table = pd.DataFrame(list_value)
# data_table

maxWind = []

for i in range(primary_count):
  new_dataFrame = data_table0[data_table0['Date'].between(starting_date[i], ending_date[i])]  #data_table0 contains 1st data
  maxWind.append(new_dataFrame['cyclonic_wind_speed'].max())

data_table['maxWind'] = maxWind
# data_table

#Msl is for mean sea level pressure
minMsl = []
maxMsl = []
meanMsl = []

for i in range(primary_count):
  minimum_pressure = np.nanmin(era_data['msl'].sel(time=slice(starting_date[i], ending_date[i])))
  maximum_pressure = np.nanmax(era_data['msl'].sel(time=slice(starting_date[i], ending_date[i])))
  mean_pressure = np.nanmean(era_data['msl'].sel(time=slice(starting_date[i], ending_date[i])))

  minMsl.append(minimum_pressure)
  maxMsl.append(maximum_pressure)
  meanMsl.append(mean_pressure)


minMsl = np.array(minMsl)
maxMsl = np.array(maxMsl)
meanMsl = np.array(meanMsl)

minMsl

data_table['minMsl(Pa)'] = minMsl
data_table['maxMsl(Pa)'] = maxMsl
data_table['meanMsl(pa)'] = meanMsl

data_table

#applying starting point rule
#if starting point's is latitude is greater than 23 then ignore
#for 1 Cyclone

first_point_remove = []
for i in range(primary_count):
  msl_hr_cyc = era_data['msl'].sel(time=slice(data_table['Starting Date'][i], data_table['Ending Date'][i]))

  cyc1_lon = []
  cyc1_lat =[]

  for j in range(len(msl_hr_cyc['time'])):
    data = msl_hr_cyc[j]
    x = np.where(data==np.nanmin(data))
    lon = data[x].longitude.data
    lat = data[x].latitude.data
    cyc1_lon.append(lon)
    cyc1_lat.append(lat)

  cyc1_lat = np.array(cyc1_lat)
  cyc1_lon = np.array(cyc1_lon)
  cyc1_lat = np.concatenate(cyc1_lat)
  cyc1_lon = np.concatenate(cyc1_lon)

  if cyc1_lat[0]>23:
    first_point_remove.append(i)


for cyc in range(len(first_point_remove)):
  data_table = data_table.drop(labels=[first_point_remove[cyc]], axis=0)
data_table = data_table.reset_index(drop=True)

print('Removed cyclone: ', first_point_remove)
data_table



# Cyclone count
second_count =len(data_table['Starting Date'])
print("Number of cyclone after applying Starting point Rule:", second_count)

#Adding duration in the dataFrame
duration = []
for i in range(second_count):
  dur = (data_table['Ending Date'][i] - data_table['Starting Date'][i])
  dur = dur / np.timedelta64(1, "h").astype("timedelta64[h]") #COnverting to hour
  duration.append(dur)
data_table['Duration'] = duration
data_table.head()

#applying cyclone duration 48 rule
dur_remove = []
for d in range(second_count):
  duration = (data_table['Ending Date'][d] - data_table['Starting Date'][d])
  duration = duration / np.timedelta64(1, 'h').astype("timedelta64[h]")

  if duration < 48:
    data_table = data_table.drop(d, axis=0)
    dur_remove.append(d)

print('Removed cyclone: ', dur_remove)
data_table = data_table.reset_index(drop=True)
data_table

# Cyclone count
third_count =len(data_table['Starting Date'])
print("Number of cyclone after applying Cyclone Duration Rule:", third_count)

print("Final Data Table")
data_table

#for all cyclone
all_lat = {}
all_lon = {}
for i in range(third_count):
  all_lat['cycLat'+str(i)]=[]
  all_lon['cycLon'+str(i)]=[]



#final
times = []
hours = []
for data in range(third_count):
  msl_hr_cyc = era_data['msl'].sel(time=slice(data_table['Starting Date'][data], data_table['Ending Date'][data]))

  for time in range(len(msl_hr_cyc['time'])):
    single_data = msl_hr_cyc[time]
    x = np.where(single_data==np.min(single_data))
    lon = single_data[x].longitude.data
    lat = single_data[x].latitude.data
    times.append(single_data.time.data)
    hours.append(time)


    all_lon['cycLon'+str(data)].append(lon)
    all_lat['cycLat'+str(data)].append(lat)

  all_lat['cycLat'+str(data)] = np.array(all_lat['cycLat'+str(data)])
  all_lon['cycLon'+str(data)] = np.array(all_lon['cycLon'+str(data)])

  all_lat['cycLat'+str(data)] = np.concatenate(all_lat['cycLat'+str(data)])
  all_lon['cycLon'+str(data)] = np.concatenate(all_lon['cycLon'+str(data)])

#applying 3.0 Degrees Rule

for y in range(third_count):
  all_lon['cycLon'+str(y)] = all_lon['cycLon'+str(y)].tolist()
  all_lat['cycLat'+str(y)] = all_lat['cycLat'+str(y)].tolist()

  n = len(all_lon['cycLon'+str(y)]) - 1
  k = 0

  while True:
    x1, y1, x2, y2 = all_lon['cycLon'+str(y)][k], all_lat['cycLat'+str(y)][k], all_lon['cycLon'+str(y)][k+1], all_lat['cycLat'+str(y)][k+1]
    distance = np.sqrt((x2-x1)**2+(y2-y1)**2)

    if distance >= 3.0:
      all_lon['cycLon'+str(y)].remove(all_lon['cycLon'+str(y)][k+1])
      all_lat['cycLat'+str(y)].remove(all_lat['cycLat'+str(y)][k+1])

      n = n-1
      if k==n:
        break
      continue
    k = k+1
    if k==n:
      break

#draw figure for all cyclone in a single figure
fig = plt.figure(figsize=(8, 10))
qq = []
#define axes
ax = plt.axes(projection=ccrs.PlateCarree())

#adding coastlines and features, skipped ocean
ax.coastlines()
ax.add_feature(cf.BORDERS, linewidth=2)
ax.add_feature(cf.LAKES)
ax.add_feature(cf.LAND)
ax.add_feature(cf.RIVERS, zorder=100)

#ploting all cyclone
for m in range(third_count):
  lon = all_lon['cycLon'+str(m)]
  qq.append(lon)
  lat = all_lat['cycLat'+str(m)]
  plt.plot(lon, lat, marker=".", markersize=4, color="teal")
  plt.scatter(lon[0], lat[0], color='r')
  plt.scatter(lon[-1], lat[-1], color='b')

#setting regional domain
ax.set_extent([80, 100, 00, 27])

#gridline for understanding lon, lat
ax.gridlines(draw_labels=True)
plt.title('Cyclone Path for cyclone',
          fontweight="bold")

lats = []

#final
for data in range(third_count):
  msl_hr_cyc = era_data['msl'].sel(time=slice(data_table['Starting Date'][data], data_table['Ending Date'][data]))

  for time in range(len(msl_hr_cyc['time'])):
    single_data = msl_hr_cyc[time]
    x = np.where(single_data==np.min(single_data))
    lat = single_data[x].latitude.data
    lats.append(lat)

nested_local = lats


flat_lat = [item for sublist in nested_local for item in sublist]

print(flat_lat)

lons = []
#final
for data in range(third_count):
  msl_hr_cyc = era_data['msl'].sel(time=slice(data_table['Starting Date'][data], data_table['Ending Date'][data]))

  for time in range(len(msl_hr_cyc['time'])):
    single_data = msl_hr_cyc[time]
    x = np.where(single_data==np.min(single_data))
    lon = single_data[x].longitude.data
    lons.append(lon)

nested_local = lons


flat_lon = [item for sublist in nested_local for item in sublist]

print(flat_lon)

x = []
for data in range(third_count):
    single_cyc = era_data['msl'].sel(time=slice(data_table['Starting Date'][data], data_table['Ending Date'][data]))
    for j in range(len(single_cyc['time'])):
        data = single_cyc[j]
        minimum_pressure0 = np.nanmin(data)
        x.append(minimum_pressure0)

o = 0
k = []
for data in range(third_count):
    single_cyc = era_data['msl'].sel(time=slice(data_table['Starting Date'][data], data_table['Ending Date'][data]))
    for j in range(len(single_cyc['time'])):
        o = data+1
        k.append(o)

x1 = []

for data in range(third_count):
    single_cyc1 = era_data['msl'].sel(time=slice(data_table['Starting Date'][data], data_table['Ending Date'][data]))
    for j in range(len(single_cyc1['time'])):
      data = single_cyc1[j]
      minimum_pressure01 = np.nanmax(data)
      x1.append(minimum_pressure01)

b = 0
h = []
for data in range(third_count):
    single_cyc = era_data['msl'].sel(time=slice(data_table['Starting Date'][data], data_table['Ending Date'][data]))
    for j in range(len(single_cyc['time'])):
        b = j+1
        h.append(b)

maxWind = []
for data in range(third_count):
    single_cyc = wind.sel(time=slice(data_table['Starting Date'][data], data_table['Ending Date'][data]))
    for j in range(len(single_cyc['time'])):
        data = single_cyc[j]
        t = np.nanmax(data)
        maxWind.append(t)



df = pd.DataFrame({'quant':k,
                   'hours':h,
                   'minMSL':x,
                   'maxwind':maxWind,
                   'lat' : flat_lat,
                   'lon' : flat_lon},

                    columns=['quant','hours','minMSL','maxwind','lat','lon'])

df

df.to_csv('data_2023.csv')