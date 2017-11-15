
## [animation of Sky SN and the ZTF survey](./two_nights_march.mp4) :
__Making the animation:__  The animation was made by running the [script](../scripts/demo_ztf.py) using [`OpSimSummary v1.2.0`](https://github.com/rbiswas4/OpSimSummary/releases/tag/v1.2.0). There are three failures at the point when the camera is in the region above 88 degrees due to problems at the poles. At that point, we did
```
cp ztf_obsHistID_032409 ztf_obsHistID_032410
cp ztf_obsHistID_032419 ztf_obsHistID_032420
cp ztf_obsHistID_032263.png ztf_obsHistID_032264.png
ffmpeg -r 15 -f image2 -s 1920x1080 -start_number 31350 -i ztf_obsHistID_%06d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p two_nights_march.mp4
```
 

