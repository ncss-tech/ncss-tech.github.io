# Survey Order, Mapping Scale, Pixel Size


1. Mapping scale (e.g. 1:24k) is somewhat less relevant than the type of mapping: order 1,2,3,4. Survey order keeps us honest and within budget. We should spend some time updating table 4-4 from https://www.nrcs.usda.gov/wps/portal/nrcs/detail/soils/ref/?cid=nrcs142p2_054254#orders to include:

a) suggested pixels sizes that may be appropriate for raster product creation, for example: 30-800m pixels for the creation of order 4 mapping products

b) suggested minimum resolvable unit for each order of mapping, this would describe filtering criteria used to remove patches of pixels (map finishing step) that are too small given the sampling density or some other constraint

A quick back-of the envelope calculation suggests that:

```
survey order | MMU (ac)  | no. 30m pixels  | no. 90m pixels
-----------------------------------------------------------
    2        | 1.5-10    |  7-45           | 1-5
    3        | 4-40      |  18-180         | 2-20
    4        | 40-60     |  180-2790       | 20-310
```

2. Pixels as component or map unit? This questions keys into two related topics: raster standards and disaggregation. SSURGO-compatible, raster-based, products should probably focus on the creation of "pixels = map unit" style maps. Single component map units would accommodate an increased level of precision--when warranted. "Pixel as component" mapping should (I think) fall into the disaggregation topic as the resulting maps won't be (immediately) SSURGO compatible.


It might be a good idea to assemble some representative landscapes from around the US and overlay couple of pieces of data:

1. SSURGO linework (various mapping orders: 2, 3, 4)
2. pixel grids (various sizes: 10m, 30m, 90m, 250m, 800m)

Then, have local experts evaluate the degree to which these grid sizes can accommodate the order of mapping used in a variety of landscapes. It could be that there will be no consensus, but I suspect that we will learn something useful. I wonder how we could do this in a controlled setting so that the experiment would not bias the results. In other works setup an experiment to determine expert judgement of pixel size / mapping order relationships.

In short, this kind of work should help us setup standards that are based on our data and our institutional knowledge... which I think is important.


## Links / References
   
   * http://spatial-analyst.net/wiki/index.php/Grid_size_calculator
   
   