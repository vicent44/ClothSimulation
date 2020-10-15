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

    void OnCollisionEnter(Collision col)
    {
        if(col.gameObject.name != "SphereLeftHand" && col.gameObject.name != "SphereRightHand" && col.gameObject.name != "New Particle")
        {
            //if(!detectedBefore)
            //{
                //this.collision = col;
                //var particle = col.gameObject.GetComponent<ParticlesBehaviour>();
                //particle.particles.isActive = false;
                //detectedBefore = true;
                //anchor = true;
                //Debug.Log("Collhand");
                //particleNum = particle.particles.I;
                //col.transform.parent = col.contacts[0].thisCollider.transform;
            //}
            _particle.Position=(col.contacts[0].point); 
            /*Vector3 d = _particle.Position - col.transform.position;
            Vector3 normalCollider = col.contacts[0].point - col.transform.position;
            Vector3.Normalize(normalCollider);

            float dot = Vector3.Dot(d,normalCollider);
            if(dot <= 0)
            {
                _particle.Position = _particle.Prev - dot * normalCollider;
            }*/
        }
    }

    void OnCollisionStay(Collision col)
    {
        if(col.gameObject.name != "SphereLeftHand" && col.gameObject.name != "SphereRightHand" && col.gameObject.name != "New Particle")
        {
            //if(anchor && detectedBefore && particleNum == col.gameObject.GetComponent<ParticlesBehaviour>().particles.I)
            //{
                //this.collision = col;
                //var particle = col.gameObject.GetComponent<ParticlesBehaviour>();
                //particle.particles.SetPosition(col.contacts[0].point);
            //}
            _particle.Position=(col.contacts[0].point); 
            /*Vector3 d = _particle.Position - col.contacts[0].point;
            Vector3 normalCollider = col.contacts[0].point - col.transform.position;
            Vector3.Normalize(normalCollider);

            float dot = Vector3.Dot(d,normalCollider);
            if(dot <= 0)
            {
                _particle.Position = _particle.Prev - dot * normalCollider;
            }*/
        }        
    }

    void OnCollisionExit(Collision col)
    {
        if(col.gameObject.name != "SphereLeftHand" && col.gameObject.name != "SphereRightHand" && col.gameObject.name != "New Particle")
        {
            //if(detectedBefore)
            //{
                //this.collision = col;
                //var particle = col.gameObject.GetComponent<ParticlesBehaviour>();
                //particle.particles.isActive = true;
                //detectedBefore = false;
            //}
        }
    }
}
