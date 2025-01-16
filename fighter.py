from vpython import *
import random
from brain import Brain, get_color

# Make scene with light blue background
scene = canvas(width=800, height=600)
scene.background = vector(0.8, 0.9, 1.0)  # Light blue using RGB values
scene.camera.pos = vector(0, 0, 100)
scene.camera.axis = vector(0, 0, -1)
scene.camera.fov = 0.1  # Smaller field of view
scene.userzoom = True   # Allow zooming
scene.userspin = True   # Allow spinning
scene.range = 50        # Set the visible range

# Replace floor line with box
box_length = 500  # -50 to 50
box_height = 12  # This will be our R value for the floor thickness
floor_y = -50  # Keep the same floor height

# Create main floor box
floor = box(pos=vector(0, floor_y - box_height/2, 0), 
           size=vector(box_length, box_height, box_height),
           color=color.green)

def random_vector(N):
    return vector(random.uniform(-N, N), random.uniform(-N, N), 0)

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
        self.radius = 12.0
        
        self.legs = []
        self.velocities = []
        for i in range(N_legs):
            leg_pos = self.pos + random_vector(9)
            leg = sphere(pos=leg_pos, radius=self.radius, color=self.color)
            self.legs.append(leg)
            self.velocities.append(vector(0, 0, 0))
            # Add to static lists
            Fighter.all_legs.append(leg)
            Fighter.all_velocities.append(self.velocities[-1])
        
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
                      radius=3.3, coils=10, color=self.color)  # Thicker springs
        return spring

    def handle_collision(self, i, j):
        # Get positions and velocities from global lists
        p1, p2 = Fighter.all_legs[i].pos, Fighter.all_legs[j].pos
        v1, v2 = Fighter.all_velocities[i], Fighter.all_velocities[j]
        
        # Calculate normal and tangential vectors
        n = (p2 - p1).norm()
        t = vector(-n.y, n.x, 0)
        
        # Project velocities
        v1n = dot(v1, n)
        v1t = dot(v1, t)
        v2n = dot(v2, n)
        v2t = dot(v2, t)
        
        # Elastic collision
        v1n_new = v2n
        v2n_new = v1n
        
        # Update velocities in global list
        Fighter.all_velocities[i] = v1n_new * n + v1t * t
        Fighter.all_velocities[j] = v2n_new * n + v2t * t
        
        # Force correct separation
        separation = 2 * self.radius
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
            # Create new brain and random color
            self.brain = Brain(2 * self.N_legs, (self.N_legs * (self.N_legs - 1)) // 2)
            self.initial_color = get_color(random.random())  # Pass random value between 0 and 1
        
        # Recreate legs and springs
        for i in range(self.N_legs):
            leg_pos = self.initial_pos + random_vector(9)
            leg = sphere(pos=leg_pos, radius=self.radius, color=self.initial_color)
            self.legs.append(leg)
            self.velocities.append(vector(0, 0, 0))
            Fighter.all_legs.append(leg)
            Fighter.all_velocities.append(self.velocities[-1])
        
        # Recreate springs
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
        
        # Get brain inputs (leg positions)
        inputs = []
        for leg in self.legs:
            inputs.extend([leg.pos.x, leg.pos.y])
        
        # Get brain outputs and update natural lengths
        outputs = self.brain.forward(inputs)
        
        # Convert outputs (0 to 1) to natural lengths (L_min to L_max)
        self.natural_lengths = [
            self.L_min + output * (self.L_max - self.L_min) 
            for output in outputs
        ]
        
        # Add gravity to forces
        forces = [vector(0, -9.8, 0) for _ in range(self.N_legs)]
        
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
            
            # Floor collision check - only if within platform bounds
            if abs(self.legs[i].pos.x) <= box_length/2:  # Use box_length instead of 50
                if (self.legs[i].pos.y < -50 + self.radius and  # Below surface
                    self.legs[i].pos.y > -50 - 5):              # But not too far below
                    self.legs[i].pos.y = -50 + self.radius  # Prevent clipping
                    self.velocities[i].y = -0.8 * self.velocities[i].y  # Elastic bounce
                    self.velocities[i].x *= 0.5  # Floor drag
                
                # Bottom collision check
                elif (self.legs[i].pos.y > -50 - box_height - self.radius and  # Above bottom surface
                      self.legs[i].pos.y < -50 - box_height + 5):             # But not too far above
                    self.legs[i].pos.y = -50 - box_height - self.radius  # Prevent clipping
                    self.velocities[i].y = -0.8 * self.velocities[i].y  # Elastic bounce
                    self.velocities[i].x *= 0.5  # Floor drag

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
    fighter1 = Fighter(vector(-20, 0, 0), 4, color.red)
    fighter2 = Fighter(vector(20, 0, 0), 3, color.blue)
    
    dt = 0.1
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
