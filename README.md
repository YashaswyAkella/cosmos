# Cosmos

A hand-built 3D map of the cosmos that runs in the browser with no libraries, no build step and no server. Open `index.html` and it just runs.

**Live:** https://YashaswyAkella.github.io/cosmos/

## What is in it

- The solar system with real orbits (JPL J2000 elements), moons, dwarf planets, the famous asteroids and comets, and spacecraft from Voyager to Webb.
- 109,000 stars, 105,000 galaxies with distances, 66,000 asteroids and comets, 16,500 satellites, 6,200 exoplanets and 3,600 nebulae and clusters, all searchable.
- Telescope views: select Hubble, Webb, Chandra or Kepler to see the patches of sky they stared at, and zoom in to their real size. The galaxies inside a patch are an impression at the real density, not the actual photograph.
- A trip planner: pick any two places and see the distance and how long the journey takes at the speed of Apollo 11, Voyager, light and more, plus the next real launch window.
- Deep time: press "Deep time ⏩" in the time bar and watch the future of the universe, with time running ten times faster every second — the planets blurring into rings, Betelgeuse exploding, Andromeda colliding with the Milky Way, the Sun becoming a red giant and then a white dwarf, the sky emptying, the last stars going out, and the black holes evaporating, out to 10^100 years. Every step is tagged as calculated, predicted or speculative.
- A size comparison: put any two things side by side to scale, from the Moon to the Milky Way.
- Real phases and rotation: every world shows its true day and night side for the simulated time, the Moon keeps its real phase, and Settings → Realistic surfaces adds Earth's continents, Jupiter's belts and the Great Red Spot, Saturn's ring divisions and shadows, the Moon's maria and Mars's dark regions. These are simplified maps drawn by hand from real geography, not photographs.
- The Milky Way itself: search for it and fly out to see the whole barred spiral from outside, with the Sun marked in the Local Arm, then keep zooming to see Andromeda beside it. No photograph of our galaxy from outside exists, so this is a model built from survey measurements of its arms and bar.
- A cosmic address for everything: every object's panel shows where it sits, from the Orion Spur through the Milky Way, the Local Group, the Virgo Supercluster and Laniakea to the edge of the observable universe, and each step is a place you can fly to.
- The deep universe: famous quasars like 3C 273 and TON 618, the Einstein Cross, the Bullet Cluster, gravitational-wave and gamma-ray-burst sources, the first galaxies, and the great walls and voids of the cosmic web, each placed by its measured distance or redshift and coloured by how long its light has travelled. The outlines of superclusters and voids are drawn at their published sizes; their real edges are irregular.
- 88 constellations by name, Sagittarius A*, pulsars, magnetars and black holes.

## Feedback wanted

This is a beta. It works best on a laptop or desktop; the phone layout is not finished yet.


## Good to know

- Positions of the planets, Moon, asteroids and comets are computed from published orbits. They are most accurate near the present: the planets stay close for a few thousand years either way, while asteroids and comets slowly drift from their true places over decades. Spacecraft and satellite positions are approximate.
- Stars, galaxies and nebulae sit at their catalogued positions and distances. Objects without a reliable distance are shown as directions on the sky.
- In Deep time, the planet orbits are real for the first thousand years; after that, and for the Sun's evolution, the Andromeda merger and everything later, the app follows published predictions and says so on each step.
- The page loads about 6 MB of catalogues on first visit, and the clock starts running at one day per second; press Now and 1× in the time bar to see the present moment.

## Controls

Drag to orbit, scroll to zoom, click to select, double-click to fly. Search at the top. The time bar at the bottom runs the clock. The horizon stays level by default; switch on Free rotation in Settings to turn in any direction, including over the poles.

## Data credits

HYG v4.0 (CC BY-SA 2.5), OpenNGC (CC BY-SA 4.0), HyperLEDA, NASA Exoplanet Archive, JPL Small-Body Database, Celestrak, Strasbourg-ESO planetary nebulae, Green's supernova remnants, Sharpless H II regions, JPL J2000 planetary elements.

Models and measurements: IAU rotation elements (Archinal et al. 2018), Meeus lunar ephemeris, Planck 2018 cosmology, Milky Way spiral arms from Reid et al. 2019, Laniakea from Tully et al. 2014, Local Group distances from McConnachie 2012, the Sun's future from Schröder & Smith 2008, and the Andromeda merger odds from Sawala et al. 2025.
