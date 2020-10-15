using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RobotCollision : MonoBehaviour
{
    private ParticlesBehaviour particle = null;
    private bool detectedBefore = false;
    private int particleNum;
    Collision collision = null;
    private bool anchor = false;

    void OnCollisionEnter(Collision col)
    {
        //if(col.gameObject.name == "New Particle")
        //{
            //if(!detectedBefore)
            //{
                this.collision = col;
                //var particle = col.gameObject.GetComponent<ParticlesBehaviour>();
                //particle.particles.isActive = false;
                //detectedBefore = true;
                //anchor = true;
                Debug.Log("Coll");
                //particleNum = particle.particles.I;
                //col.transform.parent = col.contacts[0].thisCollider.transform;
            //}
        //}
    }

    void OnCollisionStay(Collision col)
    {
        //if(col.gameObject.name == "New Particle")
       // {
            //if(anchor && detectedBefore && particleNum == col.gameObject.GetComponent<ParticlesBehaviour>().particles.I)
            //{
                this.collision = col;
                var particle = col.gameObject.GetComponent<ParticlesBehaviour>();
                //particle.particles.SetPosition(col.contacts[0].point);
            //}
        //}        
    }

    void OnCollisionExit(Collision col)
    {
        //if(col.gameObject.name == "New Particle")
        //{
            //if(detectedBefore)
            //{
                this.collision = col;
                var particle = col.gameObject.GetComponent<ParticlesBehaviour>();
                //particle.particles.isActive = true;
                //detectedBefore = false;
            //}
        //}
    }

    void Update()
    {
        /*if(Input.GetKey("p"))
        {
            Debug.Log("Unanchor the particle from the robot hand");
            anchor = false;
            //detectedBefore = false;
            collision.gameObject.GetComponent<ParticlesBehaviour>().particles.isActive = true;
            collision.gameObject.GetComponent<ParticlesBehaviour>().particles.SetPosition(collision.contacts[0].point);
            collision.contacts[0].thisCollider.transform.DetachChildren();
        }
        if(Input.GetKey("o"))
        {
            detectedBefore = false;
            Debug.Log("Now you can select another particle to move");
        }*/
    }
}
