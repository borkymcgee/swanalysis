# swanalysis

A cli tool for downloading, visualizing, and modifying full-sky lyman-a images from the SWAN instrument abord SOHO to get an estimation of solar activity.

It can:
   - Create a 3d datacube of fits data between given dates
   - Mask stars out of the dataset
   - interpolate over those masked areas
   - calculate real irradiance values
   - display an animation of the full sky image over time
     + Save that animation as a video file
     + reproject the data to a different map projection (default is mercator)

## usage:

   Run swanalysis.py -u to download the cache file. Once it completes, run swanalysis.py -h to see options and usage
   
## How does it work?
   Here i detail each of the functions swanalysis can do and how they are accomplished.
  
   ### Downloading and combining fits files:
   2d Full-sky lyman-a images for each individual day are available from the swan website in the form of fits files. Swanalysis downloads the files for all of the days over the period that you request, and combines them into one big fits file consisting of a 3d datacube of all the days stacked up.
      
   ### Masking stars
   As downloaded, the swan data is full of stars. As pretty as they are, they aren't related to the solar data we're looking for, so we cover them up with 0s. There are two phases to this step: differentiating, and floodfilling.
   The differentiation step takes advantage of the fact that the stars have quite high contrast compared with the data we want. We iterate over the data in both the x and y axes, and any pixel whose brightness is really high is assumed to be a star and gets masked out.
   The floodfill step is necessary because the differentiation step really only identifies the edges of stars. If there is a large bright region, sometimes only the edges of that region get masked, and the spot in the middle stays untouched. floodfill operates by masking any pixels that have over a certain number of neighboring pixels that are masked.
   
   ### Interpolating
   Once the stars are masked out, there are large portions of the image that are zeros not because there isn't solar activity there, but because there used to be a star (or the body of SOHO, or the sun herself) that got masked out. Interpolating serves to fill in those areas with a simple linear interpolation.
