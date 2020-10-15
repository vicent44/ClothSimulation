using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ParticlesBehaviour : MonoBehaviour
{
    public Particles particles { get {return _particle;} set {_particle = value;}}

    public bool wasKinematic { get { return _wasKinematic;} set { _wasKinematic = value;}}
    public bool moveWithMouse { get { return _moveWithMouse;} set { _moveWithMouse = value;}}

    private Particles _particle;
    private bool _wasKinematic;
    private bool _moveWithMouse;
    public Vector3 previousMousePosition = Vector3.zero;

    //Sett the inicial positions of the particles/vertex
    void Start()
    {
        transform.position = _particle.Position;
    }

    //Update the positions of the particles/vertex
    void Update()
    {
        transform.position = _particle.Position;
        
        if(Input.GetKey("p"))
        {
            transform.parent = null;
        }
    }

    //Creatio of the vertex objects
    public static ParticlesBehaviour Create(Vector3 po, float mass, int pi, int pj, GameObject baseGameObject = null)
    {
        // If a baseGameObject was passed in instantiate it, otherwise create a default
        var newGameObject = baseGameObject ? Instantiate(baseGameObject) : new GameObject();
        newGameObject.name = "New Particle";

        var newParticle = new Particles(po, mass, pi, pj);
        var newParticleBehaviour = newGameObject.AddComponent<ParticlesBehaviour>();
        var newwParticleBehaviour = newGameObject.AddComponent<Rigidbody>();
        newwParticleBehaviour.isKinematic = true;
        newwParticleBehaviour.useGravity = false;
        newParticleBehaviour._particle = newParticle;

        return newParticleBehaviour;
    }

    //Here cloth-robot collisions are handled.

    Vector3 contactPoint = Vector3.zero;

    void OnCollisionEnter(Collision col)
    {
        //First when a collision is detected, the position of the
        //contact is sabed to know the first place where the particle
        //has penetrated into the robot collider.
        if(col.gameObject.name != "SphereLeftHand" && col.gameObject.name != "SphereRightHand" && col.gameObject.name != "New Particle")
        {
            contactPoint = col.contacts[0].point;
        }
    }

    void OnCollisionStay(Collision col)
    {
        //If the collision penetrates into the collider, using the contact point 
        //the depth of the penetration is calculated and the particle
        //is pulled out of the collider by changing its position.
        if(col.gameObject.name != "SphereLeftHand" && col.gameObject.name != "SphereRightHand" && col.gameObject.name != "New Particle")
        {
            Vector3 d = _particle.Position - (contactPoint+0.01f*col.contacts[0].normal);
            float dot = Vector3.Dot(d,col.contacts[0].normal);
            if(dot <= 0)
            {
                _particle.Position = _particle.Prev - dot * col.contacts[0].normal;
            
                Vector3 normalVelocity = Vector3.Dot(col.contacts[0].normal,_particle.Velocity) * col.contacts[0].normal;
                Vector3 tangencialVelocity = _particle.Velocity - normalVelocity;
                Vector3 normalForce = Vector3.Dot(_particle.Force, col.contacts[0].normal) * col.contacts[0].normal;
                Vector3 tangencialForce = _particle.Force - normalForce;

                //_particle.Position = _particle.Position - 0.005f * (tangencialVelocity);//-frictionConstPlane*normalVelocity.magnitude*(tangencialVelocity/tangencialVelocity.magnitude)-dissipationConstPlane*normalVelocity);            
            }
        }       
    }

    void OnCollisionExit(Collision col)
    {
        //Reset of the contact point everytime that the collision finishes.
        if(col.gameObject.name != "SphereLeftHand" && col.gameObject.name != "SphereRightHand" && col.gameObject.name != "New Particle")
        {
            contactPoint = Vector3.zero;
        }
    }
}
