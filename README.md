# swanalysis

A tool for downloading, visualizing, and modifying full-sky lyman-a images from the SWAN instrument abord SOHO.

It can:
   - Create a 3d datacube of fits data between given dates
   - Mask stars out of the dataset
   - interpolate over those masked areas
   - calculate real irradiance values
   - display an animation of the full sky image over time
     + Save that animation as a video file
     + reproject the data to a different map projection (default is mercator)
