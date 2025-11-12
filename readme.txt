1) [2] Correct and automatic initialization of the particles. Every vertex in a the mesh has an associated particle, no duplicates, and in the correct coordinate system. (Hint: mesh vertices accessed via the mesh filter are in local object coordinates!).
	Done: vertices are taken from the mesh using mesh.vertices (where mesh is set to the cube mesh earlier using GetComponent<MeshFilter>().mesh) and transformed from the local object coordinates using transform.TransformPoint()

2) [3] Correct and automatic initialization of the spring configuration. No duplicates, all particles in a mesh are connected to all other particles by a single spring. The rest length is correctly computed from the initial mesh configuration.
	Done: Each spring in particles[] is connected to each other by a single spring by iterating over int j = i+1 for each particles[i]. This prevents a duplicate spring from being created as the future particle doesn't create a spring for the particles before it. The rest length is then computed by the distance between the position particle i and j.

3) [1] Correct initialization of the ground plane. Note the details of the plan come from the object representing it in the test case.
	Done: initPlane() initializes the Place variable plane by finding the object GroundPlane in the scene and using it's position (for plane.position) and rotation (for plane.normal).

4) [3] The ground contact penetration penalty springs are correctly and appropriately initialized when penetration is detected, just once during the duration of the penalty. The attach point for the penalty spring should be the nearest point on the plan from the particle at the moment the contact penetration is detection. The spring should have the correct property values and rest length. All penalty springs should have the properties k_s = 1000 and k_d = 20.
	

5) [3] The ground contact penetration penalty springs are correctly and appropriately updated during the penalty and detached when the collision is resolved.

6) [2] The vertices of the mesh are correctly updated, in the correct coordinate system, at the end of each simulation loop. Careful distinction here, you are creating "particles" mathematical objects with properties you update over time, the mesh is rendered based on where it's vertices are. You are expected to initialize a particle system using a mesh, update the particles through time, while updating the mesh vertices using the particle information.

7) [2] The particle-particle spring forces are correctly computed and the reflected force "trick" is used to reduce redundant computations of spring forces between particle pairs.

8) [2] The mesh bounds and normals are correctly updated after the mesh is modified.

9) [2] The symplectic Euler integration scheme is implemented correctly.

10) [3] The simulator loop correctly updates all particle states using the correct update callback and time. Recall the lesson on simulator loops.

11) [3] The recorded test case requires that the: blue cube has particle spring properties k_s = 200 and k_d = 0 the red cube k_s = 80 and k_d = 0.8 and the green cube k_s = 45 and k_d = 0.2.

12) You must submit a SINGLE file called <lastname>-<firstname>-a2-p2.zip (replace the <> with your information) that includes all the necessary files. You must follow the .gitignore pattern for Unity. [-2 marks if you do not follow these guidelines]
	Done

13) You must include a readme.txt file that describes in full detail which of the required elements you have completed successfully and which ones you have not. [-5 Marks if you do not include this file, and partial loss of marks for each component you do not include].
	Done: this file
