using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MouseDrag : MonoBehaviour
{


            /// <summary>
        /// The particle.
        /// </summary>
        private Particles particle;

        /// <summary>
        /// The point on the screen.
        /// </summary>
        private Vector3 screenPoint;

        /// <summary>
        /// The ray.
        /// </summary>
        private Ray ray;

        /// <summary>
        /// The hit.
        /// </summary>
        private RaycastHit hit;

        /// <summary>
        /// The late update.
        /// </summary>
        // ReSharper disable once UnusedMember.Local
        private void Update()
        {
            // Sets an anchor with right mouse click
            if (Input.GetMouseButtonDown(1))
            {
                this.ray = Camera.main.ScreenPointToRay(Input.mousePosition);
                if (Physics.Raycast(this.ray, out this.hit))
                {
                    this.hit.collider.GetComponent<Particles>().isActive = false;
                }

                this.screenPoint = Camera.main.WorldToScreenPoint(this.transform.position);
            }

            // While left mouse click is held down you can drag particles around
            if (Input.GetMouseButton(0))
                {
                    var currentScreenPoint = new Vector3(Input.mousePosition.x, Input.mousePosition.y, this.screenPoint.z);
                    var currentPosition = Camera.main.ScreenToWorldPoint(currentScreenPoint);
                    this.hit.collider.GetComponent<Particles>().Position = currentPosition;
                    this.transform.position = currentPosition;
                }

            // Unsets an anchor with middle mouse click
            if (Input.GetMouseButtonDown(2))
            {
                this.ray = Camera.main.ScreenPointToRay(Input.mousePosition);
                if (Physics.Raycast(this.ray, out this.hit))
                {
                    this.hit.collider.GetComponent<Particles>().isActive = true;
                }
            }
        }
}
