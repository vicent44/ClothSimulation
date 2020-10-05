using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RobotControl : MonoBehaviour
{
    private ParticlesBehaviour particle = null;
    private bool detectedBefore = false;
    private int particleNum;
    Collision collision = null;

    void OnCollisionEnter(Collision col)
    {
        if(col.gameObject.name == "New Particle")
        {
            if(!detectedBefore)
            {
                this.collision = col;
                var particle = col.gameObject.GetComponent<ParticlesBehaviour>();
                Debug.Log(particle.particles.isActive);
                particle.particles.isActive = false;
                detectedBefore = true;
                particleNum = particle.particles.I;

                //particle.particles.Position = col.contacts[0].point;
                //Debug.Log(col.contacts[0].thisCollider.GetComponent<SphereCollider>().radius);
                //col.contacts[0].thisCollider.transform.parent = col.transform;
                col.transform.parent = col.contacts[0].thisCollider.transform;
                //particle.transform.SetParent(col.contacts[0].thisCollider.transform);
            }
        }
    }

    void OnCollisionStay(Collision col)
    {
        if(col.gameObject.name == "New Particle")
        {
            if(detectedBefore && particleNum == col.gameObject.GetComponent<ParticlesBehaviour>().particles.I)
            {
                this.collision = col;
                var particle = col.gameObject.GetComponent<ParticlesBehaviour>();
                //particle.particles.isActive = true;
                //Debug.Log("hiStay");
                //detectedBefore = false;
                //Debug.Log(col.contacts[0].thisCollider.GetComponent<SphereCollider>().radius);
                //Vector3 distance = particle.particles.Position - transform.position;
                //float magnitudDistance = distance.magnitude;
                //Vector3.Normalize(distance);
                //Vector3 positionFinal = distance * 0.6f;
                //col.contacts[0].thisCollider.transform.parent = col.transform;
                particle.particles.Position = col.contacts[0].point;

            }
        }        
    }

    void OnCollisionExit(Collision col)
    {
        if(col.gameObject.name == "New Particle")
        {
            if(detectedBefore)
            {
                this.collision = col;
                var particle = col.gameObject.GetComponent<ParticlesBehaviour>();
                particle.particles.isActive = true;
                Debug.Log("hi");
                //detectedBefore = false;
            }
        }
    }

    void Update()
    {
        if(detectedBefore)
        {
            particle.particles.Position = transform.position;
        }
    }
}
