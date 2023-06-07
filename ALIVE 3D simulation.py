import sys
import pygame
from pygame.locals import QUIT
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import numpy as np
import random
import csv
import time as time_module

quadric = gluNewQuadric()

def create_particles(shell_count, particles_per_shell):
    particles = []
    golden_angle = np.pi * (3 - np.sqrt(5))  # Golden angle in radians

    for shell in range(1, shell_count + 1): 
        for i in range(particles_per_shell):
            phi = golden_angle * i
            y = (1 - (2 * i + 1) / (2 * particles_per_shell)) * shell
            theta = np.arccos(y / shell)  # Calculate the angle (theta) based on the y-coordinate and shell radius
            x = shell * np.cos(phi) * np.sin(theta)
            z = shell * np.sin(phi) * np.sin(theta)
            y = shell * np.cos(theta)   
            particles.append((x, y, z, shell))
    return particles

# Reference particle choice
def select_random_particle(particles):
    return random.choice(particles)

def draw_arrow(x, y, z, u, v, w, arrow_head_length, arrow_head_radius):
    arrow_length = np.sqrt(u ** 2 + v ** 2 + w ** 2)
    cylinder_length = arrow_length - arrow_head_length

    glPushMatrix()

    glTranslatef(x, y, z)

    if arrow_length != 0:
        glRotatef(np.degrees(np.arccos(w / arrow_length)), -v, u, 0)

    # Draw the cylinder part of the arrow
    gluCylinder(quadric, 0.02, 0.02, cylinder_length, 16, 1) # Arrow thickness

    # Draw the arrowhead part
    glTranslatef(0, 0, cylinder_length)
    glutSolidCone(arrow_head_radius, arrow_head_length, 16, 16)

    glPopMatrix()

def draw_particles(particles, time, shell_speeds, reference_particle, a, b, c):
    for x, y, z, shell in particles:
        speed = shell_speeds[shell - 1]  # Get the speed for the current shell

        # Calculate the unit vector components for the particle's velocity
        particle_radius = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        unit_vector_x = (x) / particle_radius
        unit_vector_y = (y) / particle_radius
        unit_vector_z = (z) / particle_radius

        # Calculate the velocity components
        u = unit_vector_x * speed * a
        v = unit_vector_y * speed * b
        w = unit_vector_z * speed * c

        # Calculate the relative velocity with respect to the reference particle
        ref_u = a * reference_particle[0] * shell_speeds[reference_particle[3] - 1] / np.sqrt(reference_particle[0] ** 2 + reference_particle[1] ** 2 + reference_particle[2] ** 2)
        ref_v = b * reference_particle[1] * shell_speeds[reference_particle[3] - 1] / np.sqrt(reference_particle[0] ** 2 + reference_particle[1] ** 2 + reference_particle[2] ** 2)
        ref_w = c * reference_particle[2] * shell_speeds[reference_particle[3] - 1] / np.sqrt(reference_particle[0] ** 2 + reference_particle[1] ** 2 + reference_particle[2] ** 2)

        rel_u = u - ref_u
        rel_v = v - ref_v
        rel_w = w - ref_w

        # Update the particle's position using its velocity components
        x_pos = u * time
        y_pos = v * time
        z_pos = w * time

        # Draw the particle
        glPushMatrix()
        glTranslatef(x_pos, y_pos, z_pos)

        # Change the color and radius of the reference particle
        if (x, y, z, shell) == reference_particle:
            glColor3f(0, 1, 1)  # Set the color to red
            glutSolidSphere(0.05, 16, 16)  # Increase the radius of partile
        else:
            glColor3f(1, 1, 1)
            glutSolidSphere(0.02, 16, 16)
        glPopMatrix()

         # Draw the velocity arrow
        if (x, y, z, shell) != reference_particle:
            arrow_head_length = 0.1
            arrow_head_radius = 0.05

            # Draw the velocity arrow only if show_arrows is True
            if show_arrows:
                glColor3f(0.75, 0, 0)  # Set the arrow color to red
                draw_arrow(
                    x_pos,
                    y_pos,
                    z_pos,
                    rel_u,
                    rel_v,
                    rel_w,
                    arrow_head_length,
                    arrow_head_radius,
                )

# Data calculation for cvs
def calculate_relative_values(particle, reference_particle, shell_speeds, time, a, b, c):
    x, y, z, shell = particle
    speed = shell_speeds[shell - 1]

    ref_x, ref_y, ref_z, ref_shell = reference_particle
    ref_speed = shell_speeds[ref_shell - 1]

    # Calculate the unit vector components for the particle's velocity
    particle_radius = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    unit_vector_x = x / particle_radius
    unit_vector_y = y / particle_radius
    unit_vector_z = z / particle_radius

    # Calculate the unit vector components for the reference particle's velocity
    ref_particle_radius = np.sqrt(ref_x ** 2 + ref_y ** 2 + ref_z ** 2)
    ref_unit_vector_x = ref_x * a / ref_particle_radius
    ref_unit_vector_y = ref_y * b / ref_particle_radius
    ref_unit_vector_z = ref_z * c / ref_particle_radius

    # Calculate velocity components
    u = unit_vector_x * speed * a
    v = unit_vector_y * speed * b
    w = unit_vector_z * speed * c

    ref_u = ref_unit_vector_x * ref_speed
    ref_v = ref_unit_vector_y * ref_speed
    ref_w = ref_unit_vector_z * ref_speed

    # Calculate current positions
    x_curr = u * time
    y_curr = v * time
    z_curr = w * time

    ref_x_curr = ref_u * time
    ref_y_curr = ref_v * time
    ref_z_curr = ref_w * time

    # Calculate relative position components
    rel_x = x_curr - ref_x_curr
    rel_y = y_curr - ref_y_curr
    rel_z = z_curr - ref_z_curr

    # Calculate relative distance
    rel_dist = np.sqrt(rel_x**2 + rel_y**2 + rel_z**2)

    # Calculate relative velocity components
    rel_u = u - ref_u
    rel_v = v - ref_v
    rel_w = w - ref_w

    # Calculate relative velocity
    rel_vel = (np.sqrt(rel_u**2 + rel_v**2 + rel_w**2)) * 1000  # Scaled for clarity

    # Calculate relative Hubble value
    rel_hubble = rel_vel / rel_dist

    return rel_dist, rel_vel, rel_hubble

def save_data_to_csv(data, filename):
    with open(filename, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Time (minutes)","Shell #", "Particle #", "Relative Distance", "Relative Velocity", "Relative Hubble Value"])
        # Increase the precision for floating-point values
        np.set_printoptions(precision=16, suppress=True)
        csv_writer.writerows(data)

show_arrows = False 

# Main Loop
def main():
    pygame.init()
    display = (1500, 900)
    pygame.display.set_mode(display, pygame.DOUBLEBUF | pygame.OPENGL)

    particles = create_particles(5,100) # Shells, particle per shells <------------------------

    elapsed_time = 0
    data = []
    start_time = time_module.time()

    global show_arrows

    # Select a random particle as the reference frame
    reference_particle = select_random_particle(particles)

    clock = pygame.time.Clock()
    time = 0.0

    # Initialize rotation angles
    x_rotation = 0
    y_rotation = 0

    # Initialize zoom level
    zoom = -20

    # Initialize mouse button and pause states
    mouse_button_down = False
    paused = False

    while True:
        # Capture the time since the last data recording
        current_time = time_module.time()
        if current_time - start_time >= 1:  #Adjust time interval that data is collected, 1 = 1 second
            start_time = current_time
            elapsed_time += 1
            for i, particle in enumerate(particles):
                rel_dist, rel_vel, rel_hubble = calculate_relative_values(particle, reference_particle, shell_speeds, time, a, b, c)
                data.append((elapsed_time, particle[3], i, rel_dist, rel_vel, rel_hubble))
        for event in pygame.event.get():
            if event.type == QUIT:
                save_data_to_csv(data, "particle_data.csv")
                pygame.quit()
                sys.exit()
            # Update mouse button state
            elif event.type == pygame.MOUSEBUTTONDOWN:
                if event.button == 1:  # Left mouse button
                    mouse_button_down = True
            elif event.type == pygame.MOUSEBUTTONUP:
                if event.button == 1:  # Left mouse button
                    mouse_button_down = False
            # Capture mouse motion events
            elif event.type == pygame.MOUSEMOTION and mouse_button_down:
                x_rel, y_rel = event.rel
                x_rotation += y_rel * 0.1
                y_rotation += x_rel * 0.1
            # Toggle pause state with spacebar
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE:
                    paused = not paused
                # Change the reference particle with "R" key
                elif event.key == pygame.K_r:
                    reference_particle = select_random_particle(particles)
                # Zoom in with + key
                elif event.key in (pygame.K_EQUALS, pygame.K_KP_EQUALS):
                    zoom += 0.5
                # Zoom out with - key
                elif event.key in (pygame.K_MINUS, pygame.K_KP_MINUS):
                    zoom -= 0.5
                # Toggle arrows visibility with "A" key
                elif event.key == pygame.K_a: 
                    show_arrows = not show_arrows
        current_time = pygame.time.get_ticks()

        glLoadIdentity()
        gluPerspective(45, (display[0] / display[1]), 0.1, 50.0)
        glTranslatef(0.0, 0.0, zoom)
        glRotatef(x_rotation, 1, 0, 0)
        glRotatef(y_rotation, 0, 1, 0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        shell_speeds = range(1,101) # Shell speed range
        a,b,c = [0.1,0.000001,0.1] # Use this to change the overall speeds

        draw_particles(particles, time, shell_speeds, reference_particle, a, b, c)

        if not paused:
            time += 0.01
        pygame.display.flip() 
        clock.tick(60)

if __name__ == '__main__':
    main()
