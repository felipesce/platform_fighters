from vpython import *
import random
from brain import Brain, get_color

dt = 0.06
N1, N2 = 3, 3
P1 = vec(-25,0,0)
P2 = -P1

# Make scene with light blue background
scene = canvas(width=800, height=600)
scene.background = vector(0.8, 0.9, 1.0)  # Light blue using RGB values
scene.camera.pos = vector(0, 0, 100)
scene.camera.axis = vector(0, 0, -1)
scene.userzoom = True   # Allow zooming
scene.userspin = True   # Allow spinning
scene.range = 250        # Set the visible range

# Define platform dimensions
box_length = 500  # Define the base size
platform_radius = box_length/2  # Using same scale as before
floor_y = -50  # Keep the same floor height
box_height = 12  # Keep same thickness

# Create main floor as a cylinder
floor = cylinder(pos=vector(0, floor_y - box_height/2, 0),
                axis=vector(0, box_height, 0),
                radius=platform_radius,
                color=color.green)

def random_vector(N):
    # Modified to include z-coordinate but ensure CM stays at z=0
    return vector(random.uniform(-N, N), random.uniform(-N, N), random.uniform(-N, N))

class Fighter:
    all_legs = []  # Static list to store all legs
    all_velocities = []  # Static list to store all velocities

    def __init__(self, pos, N_legs, color=color.red):
        self.initial_pos = pos  # Store initial position for reset
        self.initial_color = color  # Store initial color for reset
        self.pos = pos
        self.N_legs = N_legs
        self.color = color
        self.k = 2.0
        self.L = 67.0
        self.damping = 1.0
        self.radius = 20.0
        
        self.legs = []
        self.velocities = []
        # Modified leg initialization to maintain z=0 center of mass
        total_z = 0
        for i in range(N_legs):
            leg_pos = self.pos + random_vector(20)
            total_z += leg_pos.z
            leg = sphere(pos=leg_pos, radius=self.radius, color=self.color)
            self.legs.append(leg)
            self.velocities.append(vector(0, 0, 0))
            
        # Adjust z positions to center CM at z=0
        z_offset = total_z / N_legs
        for leg in self.legs:
            leg.pos.z -= z_offset
            # Add to static lists
            Fighter.all_legs.append(leg)
            Fighter.all_velocities.append(vector(0, 0, 0))
        
        self.springs = []
        for i in range(N_legs):
            for j in range(i+1, N_legs):
                spring = self.connect_legs(i, j)
                self.springs.append((i, j, spring))

        # Add L_min and L_max as class attributes
        self.L_min = 60.0
        self.L_max = 120.0
        self.L = (self.L_max + self.L_min) / 2  # Set default L to middle value
        
        # Calculate number of inputs: 2 for each leg's position (x,y)
        n_inputs = 2 * N_legs
        # Number of outputs = number of spring connections
        n_outputs = (N_legs * (N_legs - 1)) // 2
        
        self.brain = Brain(n_inputs, n_outputs)
        self.natural_lengths = [self.L] * n_outputs  # Initialize all springs with default length

    def connect_legs(self, i, j):
        spring = helix(pos=self.legs[i].pos, axis=self.legs[j].pos - self.legs[i].pos,
                      radius=10.0, coils=10, color=self.color)  # Thicker springs
        return spring

    def handle_collision(self, i, j):
        # Modified for full 3D collision handling
        p1, p2 = Fighter.all_legs[i].pos, Fighter.all_legs[j].pos
        v1, v2 = Fighter.all_velocities[i], Fighter.all_velocities[j]
        
        # Calculate normal vector in 3D
        n = (p2 - p1).norm()
        
        # Calculate two perpendicular tangent vectors
        # First tangent vector
        if abs(n.x) > abs(n.y):
            t1 = vector(-n.z, 0, n.x).norm()
        else:
            t1 = vector(0, -n.z, n.y).norm()
        # Second tangent vector (cross product)
        t2 = cross(n, t1)
        
        # Project velocities
        v1n = dot(v1, n)
        v1t1 = dot(v1, t1)
        v1t2 = dot(v1, t2)
        v2n = dot(v2, n)
        v2t1 = dot(v2, t1)
        v2t2 = dot(v2, t2)
        
        # Elastic collision calculations
        e = 0.8
        m1 = m2 = 1.0
        
        v1n_new = (v1n * (m1 - e*m2) + v2n * m2 * (1 + e)) / (m1 + m2)
        v2n_new = (v2n * (m2 - e*m1) + v1n * m1 * (1 + e)) / (m1 + m2)
        
        # Update velocities with all components
        Fighter.all_velocities[i] = v1n_new * n + v1t1 * t1 + v1t2 * t2
        Fighter.all_velocities[j] = v2n_new * n + v2t1 * t1 + v2t2 * t2
        
        # Update positions to prevent sticking
        separation = 2.1 * self.radius
        midpoint = (p1 + p2) / 2
        Fighter.all_legs[i].pos = midpoint - n * separation/2
        Fighter.all_legs[j].pos = midpoint + n * separation/2

    def reset(self, new_brain=False):
        """Reset fighter position and optionally create new brain
        new_brain: if True, create new brain and random color; if False, just reposition"""
        
        # Remove old legs and springs from visualization
        for leg in self.legs:
            leg.visible = False
            Fighter.all_legs.remove(leg)
        for _, _, spring in self.springs:
            spring.visible = False
        
        # Clear instance lists
        self.legs.clear()
        self.velocities.clear()
        self.springs.clear()
        
        if new_brain:
            # Create new brain and select from fixed colors
            self.brain = Brain(2 * self.N_legs, (self.N_legs * (self.N_legs - 1)) // 2)
            pink = vector(1.0, 0.75, 0.8)  # Define pink using RGB values
            colors = [color.red, color.green, color.blue, color.cyan, 
                     color.magenta, color.yellow, pink]
            # Remove current color of other fighters from available colors
            available_colors = [c for c in colors if not any(
                isinstance(f, Fighter) and f != self and 
                abs(c.x - f.initial_color.x) < 0.01 and
                abs(c.y - f.initial_color.y) < 0.01 and
                abs(c.z - f.initial_color.z) < 0.01
                for f in globals().values()
            )]
            self.initial_color = random.choice(available_colors)
        
        # Recreate legs and springs with current color
        for i in range(self.N_legs):
            leg_pos = self.initial_pos + random_vector(9)
            leg = sphere(pos=leg_pos, radius=self.radius, color=self.initial_color)
            self.legs.append(leg)
            self.velocities.append(vector(0, 0, 0))
            Fighter.all_legs.append(leg)
            Fighter.all_velocities.append(self.velocities[-1])
        
        # Recreate springs with same color as legs
        for i in range(self.N_legs):
            for j in range(i+1, self.N_legs):
                spring = self.connect_legs(i, j)
                self.springs.append((i, j, spring))

    def update(self, dt):
        # Add death check at the start of update
        for leg in self.legs:
            if leg.pos.y < -250:
                self.reset()
                return
        
        # Get inputs from leg positions (x,y coordinates)
        inputs = []
        for leg in self.legs:
            inputs.extend([leg.pos.x, leg.pos.y])
        
        # Get outputs from brain using forward() method
        outputs = self.brain.forward(inputs)
        
        # Convert outputs (0 to 1) to natural lengths (L_min to L_max)
        self.natural_lengths = [
            self.L_min + output * (self.L_max - self.L_min) 
            for output in outputs
        ]
        
        # Add gravity and air drag to forces
        gamma = 0.0043  # Air drag coefficient
        forces = []
        for i in range(self.N_legs):
            gravity = vector(0, -9.8, 0)
            air_drag = -gamma * self.velocities[i]
            forces.append(gravity + air_drag)

        # Calculate spring forces (updated to use individual natural lengths)
        for idx, (i, j, spring) in enumerate(self.springs):
            displacement = self.legs[j].pos - self.legs[i].pos
            current_length = mag(displacement)
            direction = displacement.norm()
            
            # Use specific natural length for this spring
            L = self.natural_lengths[idx]
            force = self.k * (current_length - L) * direction
            
            forces[i] += force
            forces[j] -= force
            
            spring.pos = self.legs[i].pos
            spring.axis = displacement

        # Update velocities and positions
        for i in range(self.N_legs):
            self.velocities[i] += forces[i] * dt
            self.velocities[i] *= self.damping
            self.legs[i].pos += self.velocities[i] * dt
            
            # Floor collision check - only if within platform radius
            radial_dist = sqrt(self.legs[i].pos.x**2 + self.legs[i].pos.z**2)
            if radial_dist <= platform_radius:  # Check if within circle
                if (self.legs[i].pos.y < -50 + self.radius and  # Below surface
                    self.legs[i].pos.y > -50 - 5):              # But not too far below
                    self.legs[i].pos.y = -50 + self.radius  # Prevent clipping
                    self.velocities[i].y = -0.8 * self.velocities[i].y  # Elastic bounce
                    self.velocities[i].x *= 0.5  # Floor drag
                    self.velocities[i].z *= 0.5  # Floor drag in z direction too
                
                # Bottom collision check
                elif (self.legs[i].pos.y > -50 - box_height - self.radius and  # Above bottom surface
                      self.legs[i].pos.y < -50 - box_height + 5):             # But not too far above
                    self.legs[i].pos.y = -50 - box_height - self.radius  # Prevent clipping
                    self.velocities[i].y = -0.8 * self.velocities[i].y  # Elastic bounce
                    self.velocities[i].x *= 0.5  # Floor drag
                    self.velocities[i].z *= 0.5  # Floor drag in z direction too

        # Check collisions with ALL legs
        for i, leg in enumerate(self.legs):
            my_leg_index = Fighter.all_legs.index(leg)
            # Check against all other legs
            for j, other_leg in enumerate(Fighter.all_legs):
                if other_leg != leg:  # Don't check against self
                    p1, p2 = leg.pos, other_leg.pos
                    distance = mag(p2 - p1)
                    if distance < 2 * self.radius:
                        # Handle collision using existing method but with global indices
                        self.handle_collision(my_leg_index, j)

    
if __name__ == "__main__":
    # Create two fighters further apart
    fighter1 = Fighter(P1, N1, color.red)
    fighter2 = Fighter(P2, N2, color.blue)
    
    while True:
        rate(30)
        
        # Store initial states before updates
        fighter1_died = any(leg.pos.y < -250 for leg in fighter1.legs)
        fighter2_died = any(leg.pos.y < -250 for leg in fighter2.legs)
        
        # If either dies, reset both but only regenerate the dead one
        if fighter1_died or fighter2_died:
            fighter1.reset(new_brain=fighter1_died)
            fighter2.reset(new_brain=fighter2_died)
        else:
            fighter1.update(dt)
            fighter2.update(dt)
