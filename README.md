#Inspiration
The Sensel Morph's ability to sense varying pressure in a 2-D plane and its exceptional multitouch capability naturally lends it to interacting with spatial simulations. We decided that using the multi-touch capabilities and pressure sensitivity for a physics simulation would be interesting.

#What it does
Our simulation allows for intuitive and simultaneous creation of massive bodies using pressure for size and direction of dragging for velocity. It simulates attractive interactions as well as collisions. With it, the user can create a simple binary system or have fun watching tons of spheres collide around the screen at different angles.

#How we built it
We used the Processing language and the Sensel API to create this program.

#Challenges we ran into
We had some issues with learning Processing for the first time, as well as some issues with interfacing with the Sensel API. We also spent a lot of time figured out how to calculate 3D collisions. In addition, drawing the spheres with proper lighting was a challenge, with very little documentation online to refer to.

#Accomplishments that we're proud of
We learned how to interface with Sensel Morph, as well as how to use Processing. We also wrote complicated algorithms to calculate 3D elastic collisions and rotate points in 3D space.

#What I learned
Abstract is not always better. It may have been more expedient to write this in c++ and OpenGL, as that would have allowed us to avoid some of the known issues inherent in Processing (although many of these have been fixed in 3.0+). However, learning to create a program through Processing, which was written for visual programs, was nonetheless valuable and beneficial to our continued education in computer science.

#What's next for Sensel Morph Planets
While there are many future possibilities for this project, the immediately succeeding ones include the following: We could add the ability to move around and scale the system, rotate the system, and add planets in different areas in that manner. We also may add varying elasticity for planets (all currently have uniform elasticity). We could also add curves showing the predicted and past paths of planets. We could also add a more user-friendly interface to toggle all of these options.
