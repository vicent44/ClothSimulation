using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MoveTarget : MonoBehaviour
{
    public Transform target;
    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        if(Input.GetKey("left"))
        {
            target.position +=(new Vector3(0.01f, 0f, 0f));
        }
        if(Input.GetKey("right"))
        {
            target.position +=(new Vector3(-0.01f, 0f, 0f));
        }
        if(Input.GetKey("up"))
        {
            target.position +=(new Vector3(0f, 0f, 0.01f));
        }
        if(Input.GetKey("down"))
        {
            target.position +=(new Vector3(0f, 0f, -0.01f));
        }
        if(Input.GetKey("w"))
        {
            target.position +=(new Vector3(0f, 0.01f, 0f));
        }
        if(Input.GetKey("s"))
        {
            target.position +=(new Vector3(0f, -0.01f, 0f));
        }        
    }
}
