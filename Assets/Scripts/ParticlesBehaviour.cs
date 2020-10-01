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
}
