using System.Collections;
using System.Collections.Generic;
using UnityEngine;

// Check this out we can require components be on a game object!
[RequireComponent(typeof(MeshFilter))]

public class BParticleSimMesh : MonoBehaviour
{
    public struct BSpring
    {
        public float kd;                        // damping coefficient
        public float ks;                        // spring coefficient
        public float restLength;                // rest length of this spring
        public int attachedParticle;            // index of the attached other particle (use me wisely to avoid doubling springs and sprign calculations)
    }

    public struct BContactSpring
    {
        public float kd;                        // damping coefficient
        public float ks;                        // spring coefficient
        public float restLength;                // rest length of this spring (think about this ... may not even be needed o_0)
        public Vector3 attachPoint;             // the attached point on the contact surface
    }

    public struct BParticle
    {
        public Vector3 position;                // position information
        public Vector3 velocity;                // velocity information
        public float mass;                      // mass information
        public BContactSpring contactSpring;    // Special spring for contact forces
        public bool attachedToContact;          // is this particle currently attached to a contact (ground plane contact)
        public List<BSpring> attachedSprings;   // all attached springs, as a list in case we want to modify later fast
        public Vector3 currentForces;           // accumulate forces here on each step        
    }

    public struct BPlane
    {
        public Vector3 position;                // plane position
        public Vector3 normal;                  // plane normal
    }

    public float contactSpringKS = 1000.0f;     // contact spring coefficient with default 1000
    public float contactSpringKD = 20.0f;       // contact spring daming coefficient with default 20

    public float defaultSpringKS = 100.0f;      // default spring coefficient with default 100
    public float defaultSpringKD = 1.0f;        // default spring daming coefficient with default 1

    public bool debugRender = false;            // To render or not to render


    /*** 
     * I've given you all of the above to get you started
     * Here you need to publicly provide the:
     * - the ground plane transform (Transform)
     * - handlePlaneCollisions flag (bool)
     * - particle mass (float) 
     * - useGravity flag (bool)
     * - gravity value (Vector3)
     * Here you need to privately provide the:
     * - Mesh (Mesh)
     * - array of particles (BParticle[])
     * - the plane (BPlane)
     ***/

    // public variables
    // public Transform groundPlaneTransform; ???
    public bool handlePlaneCollisions = true;
    public const float particle_mass = 1.0f; // kg 
    public bool useGravity = true;
    public Vector3 gravity = new Vector3(0.0f*particle_mass, -9.8f*particle_mass, 0.0f*particle_mass);

    // private variables
    private Mesh mesh;
    private const int n = 1; // size of array
    private BParticle[] particles = new BParticle[n]; 
    private BPlane plane;
    

    private Vector3[] new_positions = new Vector3[n]; // x_{i, new}
    private Vector3[] new_velocitys = new Vector3[n]; // v_{i, new}



    /// <summary>
    /// Init everything
    /// HINT: in particular you should probbaly handle the mesh, init all the particles, and the ground plane
    /// HINT 2: I'd for organization sake put the init particles and plane stuff in respective functions
    /// HINT 3: Note that mesh vertices when accessed from the mesh filter are in local coordinates.
    ///         This script will be on the object with the mesh filter, so you can use the functions
    ///         transform.TransformPoint and transform.InverseTransformPoint accordingly 
    ///         (you need to operate on world coordinates, and render in local)
    /// HINT 4: the idea here is to make a mathematical particle object for each vertex in the mesh, then connect
    ///         each particle to every other particle. Be careful not to double your springs! There is a simple
    ///         inner loop approach you can do such that you attached exactly one spring to each particle pair
    ///         on initialization. Then when updating you need to remember a particular trick about the spring forces
    ///         generated between particles. 
    /// </summary>
    void Start(){
        initPlane();
        initParticles();
    }


    private void initPlane(){
        plane = new BPlane();
        plane.position = new Vector3(0.0f, -3.36f, 0.0f);
        plane.normal = new Vector3(0.0f, 1.0f, 0.0f);
    }

    private void initParticles(){
        BContactSpring cS= new BContactSpring();
        cS.kd = contactSpringKD; 
        cS.ks = contactSpringKS;
        cS.restLength = 0.0f;
        cS.attachPoint = plane.position;


        for(int i = 0; i < n; i++){
            BParticle p = new BParticle();
            /*** create new BParticle
                public Vector3 position;                // position information
                public Vector3 velocity;                // velocity information
                public float mass;                      // mass information
                public BContactSpring contactSpring;    // Special spring for contact forces
                public bool attachedToContact;          // is this particle currently attached to a contact (ground plane contact)
                public List<BSpring> attachedSprings;   // all attached springs, as a list in case we want to modify later fast
                public Vector3 currentForces; 
            ***/
            //p.position = [position from Mesh *changing coordinate system*]
            p.velocity = new Vector3(0.0f, 0.0f, 0.0f);
            p.mass = particle_mass;
            p.contactSpring = cS;
            //p.attachedToContact = t/f;
            particles[i].currentForces = (useGravity)? gravity : new Vector3(0.0f, 0.0f, 0.0f);
            particles[i] = p;
        }
    }

    /*** BIG HINT: My solution code has as least the following functions
     * InitParticles()
     * InitPlane()
     * UpdateMesh() (remember the hint above regarding global and local coords)
     * ResetParticleForces()
     * ...
     ***/


    public void fixedUpdate(){
        
        // First we calculate the Forces acting on each particle
        resetParticleForces();

        // Next, we use that forces to calculate all the particles new velocities and positions using the Symplectic Euler integration scheme (slide 45)
        // *We don't update the actual particle information until they have been calculated for all the particles*
        for(int i = 0; i < n; i++) {
            new_velocitys[i] = particles[i].velocity + Time.fixedDeltaTime*(particles[].currentForces/particles[i].mass); // v_{i, new} = v_i + dt*(F/m)
            new_positions[i] = particles[i].position + Time.fixedDeltaTime*new_velocitys[i]; // x_{i, new} = x_i + dt*v_{i, new}
        }

        // Once we are done calculating each particles new velocity and position, we update them
        for(int i = 0; i < n; i++){
            particles[i].velocity = new_velocitys[i]; 
            particles[i].position = new_positions[i]; 
        }

        // Finally, we update the mesh with the new positions
        updateMesh();
    }

    private void updateMesh(){
        for(int i = 0; i < n; i++){
            // don't forget to change back into Mesh Coordinate system!
        }
    }

    private void resetParticleForces(){
        /***
            There are three types of forces that are acting on these particles:
                1) gravity
                2) ground penetration penalty forces
                3) particle-particle spring forces
        ***/

        for(int i = 0; i < n; i++){
            BParticle curr_particle = particles[i];

            // 1) gravity *if* useGravity is toggled on
            curr_particle.currentForces = (useGravity)? gravity : new Vector3(0.0f, 0.0f, 0.0f);

            // 2) ground penetration penalty forces
                // -k_s((x_p - x_g)dot n)n - k_d*v_p
                // we only want to do this one if a collision is detected
                // how we detect a collision?
            
            // first, we update the attach point to be directly under the position of the current particle (using the y-coor of the plane)
            curr_particle.contactSpring.attachPoint = new Vector3(curr_particle.position.x, plane.position.y, curr_particle.position.z)
            // Next, we detect if there is a collision using the dot product of the normal vector of the plane and the vector from the position of the particle to the attach point
            Vector3 d = Vector3.Dot(curr_particle.position - curr_particle.contactSpring.attachPoint, plane.normal) // (x_p-x_g)dot n
            if(d < 0.0){ // (x_p-x_g)dot n < 0
                // if there is a collision, we calculate the ground contact penetration penalty spring equation given in the assignment
                BContactSpring curr_contact_spring = curr_particle.contactSpring;
                Vector3 ground_penalty = -1*curr_contact_spring.k_s*d*plane.normal - curr_contact_spring.k_d*curr_particle.velocity;// -k_s((x_p - x_g)dot n)n - k_d*v_p
                // and add it to the current forces variable of our current particle
                curr_particle.currentForces += ground_penalty;
            }


            // add any particle-particle spring forces that have not already been calculated
            // 3) particle-particle spring forces
            // we only calculate the forces that have not been calculated previously
            for(int j = 0; j < particles[i].attachedSprings.Count; j++){
                BSpring curr_spring = particles[i].attachedSprings[j];

                if (curr_spring.attachedParticle > i) {
                    // if the index of the other particle in the string is less than the index of our current particle, we have calculated the spring force already
                    // if not, we calculate f_ij and add it to force for both the current and corresponding particles
                    Vector3 f_ij = springForce(curr_particle, curr_spring); 
                    curr_particle.currentForces += f_ij;
                    particles[curr_spring.attachedParticle].currentForces -= f_ij;
                }

            }
        }
    }

    private Vector3 springForce(BParticle i, BSpring s){
        // this is just a helper function for readability, it is only called per particle-particle spring
        Vector3 f = new Vector3(0.0f, 0.0f, 0.0f);
        BParticle attached = particles[s.attachedParticle];
        Vector3 differenceVector = i.position - attached.position;
        Vector3 dVnormalized = differenceVector/differenceVector.magnitude;
        f += s.ks*(s.restLength - differenceVector.magnitude)*dVnormalized; 
        f += s.kd*(Vector3.Dot(i.velocity - attached.velocity, dVnormalized)*dVnormalized);
        return f;
    }

    /// <summary>
    /// Draw a frame with some helper debug render code
    /// </summary>
    public void Update()
    {
        /* This will work if you have a correctly made particles array
        if (debugRender)
        {
            int particleCount = particles.Length;
            for (int i = 0; i < particleCount; i++)
            {
                Debug.DrawLine(particles[i].position, particles[i].position + particles[i].currentForces, Color.blue);

                int springCount = particles[i].attachedSprings.Count;
                for (int j = 0; j < springCount; j++)
                {
                    Debug.DrawLine(particles[i].position, particles[particles[i].attachedSprings[j].attachedParticle].position, Color.red);
                }
            }
        }
        */
    }
}
