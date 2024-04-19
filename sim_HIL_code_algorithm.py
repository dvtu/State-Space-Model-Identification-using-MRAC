from machine import ADC, Pin, PWM
import ulab
from ulab import numpy as np
from machine import Timer
import uio

x_m_prev = np.array([0, 0])
A_m_prev = np.array([[0, 0], [0, 0]])
B_m_prev = np.array([0, 0])

I_max = 0.5
theta_max = 5
u_max = 12
voltage = 1
pwm = PWM(Pin(22))
pwm.freq(1000)
duty = int((voltage * 65535) / u_max)
pwm.duty_u16(duty)
adc_ch1 = ADC(Pin(27))
adc_ch2 = ADC(Pin(28))
read_pin = Pin(17, Pin.IN)
led = Pin(25, Pin.OUT)

# Create a memory stream
stream = uio.StringIO()
# Initialize sample count
sample_count = 0
sample_max = 1000
def algorithm(pin):
    global x_m_prev, A_m_prev, B_m_prev, sample_count
    
    current = (adc_ch2.read_u16() / 65535) * I_max
    speed = (adc_ch1.read_u16() / 65535) * theta_max
    
    dt = 0.01    
    P = np.array([[1, 1],[1, 3]])
    Gamma_a = np.array([[100, 0],[0, 100]])
    Gamma_b = np.array([[100, 0],[0, 100]])
    
    x = np.array([current, speed])    
 
    x_m = x_m_prev + (np.dot(A_m_prev, x_m_prev) + B_m_prev * voltage) * dt
    
    e = x - x_m
    
    term1 = np.dot(np.dot(Gamma_a, P), e.reshape((2, 1))) #2d-array 2x1
    term2 = np.dot(term1, np.array([x_m]))
    A_m = A_m_prev + term2 * dt
    term3 = np.dot(np.dot(Gamma_b, P), e.reshape((2, 1))) #2d-array 2x1
    B_m = B_m_prev + term3.flatten() * voltage * dt
    
    x_m_prev = x_m
    A_m_prev = A_m
    B_m_prev = B_m
    
    # Append data to the memory stream
    if sample_count < sample_max:
        stream.write("{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(A_m[0,0], A_m[0,1], A_m[1,0], A_m[1,1], B_m[0], B_m[1]))     
    sample_count += 1  # Increment sample count
    if sample_count >= sample_max:
        with open('data.txt', 'a') as f:
            f.write(stream.getvalue())
        stream.seek(0)  # Reset the stream position     
    print(sample_count)
# Configure interrupt on Pin 17 for rising edge
read_pin.irq(trigger=Pin.IRQ_RISING, handler=algorithm)
while True:
    pass
