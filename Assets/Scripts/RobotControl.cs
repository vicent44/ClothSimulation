using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RobotControl : MonoBehaviour
{
    private ParticlesBehaviour particle;
    private bool detectedBefore = false;
    private int particleNum;

    void OnCollisionEnter(Collision col)
    {
        if(col.gameObject.name == "New Particle")
        {
            if(!detectedBefore)
            {
                //Debug.Log("col.transform");
                var particle = col.gameObject.GetComponent<ParticlesBehaviour>();
                Debug.Log(particle.particles.isActive);
                particle.particles.isActive = false;
                //Debug.Log(particle.GetType());
                //Debug.Log(particle.particles.isActive);
                detectedBefore = true;
                particleNum = particle.particles.I;
                //oldPos = 
                //FixedJoint joint = gameObject.AddComponent<FixedJoint>();
                //joint.anchor = col.contacts[1].point;
                //joint.connectedBody = col.contacts[1].otherCollider.transform.GetComponentInParent<Rigidbody>();
                //joint.enableCollision = false;
            }
        }
    }

    void OnCollisionStay(Collision col)
    {
        if(col.gameObject.name == "New Particle")
        {
            if(detectedBefore && particleNum == col.gameObject.GetComponent<ParticlesBehaviour>().particles.I)
            {
                var particle = col.gameObject.GetComponent<ParticlesBehaviour>();
                //particle.particles.isActive = true;
                //Debug.Log("hiStay");
                //detectedBefore = false;
                Debug.Log(col.contacts[1].thisCollider.GetComponent<SphereCollider>().radius);
                //Vector3 distance = particle.particles.Position - transform.position;
                //float magnitudDistance = distance.magnitude;
                //Vector3.Normalize(distance);
                //Vector3 positionFinal = distance * 0.6f;
                particle.particles.Position = col.contacts[0].point;//transform.position + positionFinal;

            }
        }        
    }

    void OnCollisionExit(Collision col)
    {
        if(col.gameObject.name == "New Particle")
        {
            if(detectedBefore)
            {
                var particle = col.gameObject.GetComponent<ParticlesBehaviour>();
                particle.particles.isActive = true;
                Debug.Log("hi");
                detectedBefore = false;
            }
        }
        //particle.particles.isActive = true;
    }
    // Start is called before the first frame update
    /*void Start()
    {
        if(particle != null)
        {
            particle.particles.isActive = false;
            particle.particles.Position = transform.position;
        }
    }

    // Update is called once per frame
    void Update()
    {
        if(particle != null)
        {
            particle.particles.isActive = false;
            particle.particles.Position = transform.position;
        }
    }*/
}
