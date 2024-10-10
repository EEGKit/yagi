use crate::modem::modem::*;

#[derive(Debug, Clone)]
pub(super) struct Pi4Dqpsk {
    theta: f32,    // phase state
}

impl Pi4Dqpsk {
    pub(super) fn reset(&mut self) {
        self.theta = 0.0;
    }
}

impl Modem {
    pub(super) fn new_pi4dqpsk() -> Result<Self> {
        let mut modem = Modem::_new(2, ModulationScheme::Pi4Dqpsk)?;
        let data = Pi4Dqpsk { theta: 0.0 };
        modem.data = Some(ModemData::Pi4Dqpsk(data));
        modem.demodulate_soft_func = Some(Self::demodulate_soft_pi4dqpsk);
        Ok(modem)
    }

    pub(super) fn modulate_pi4dqpsk(&mut self, sym_in: u32) -> Result<Complex32> {
        if let ModemData::Pi4Dqpsk(pi4dqpsk) = self.data.as_mut().unwrap() {
            let d_theta = match sym_in {
                0 => 0.25 * PI,
                1 => 0.75 * PI,
                2 => -0.25 * PI,
                3 => -0.75 * PI,
                _ => return Err(Error::Config("invalid input symbol".to_string())),
            };

            pi4dqpsk.theta += d_theta;

            if pi4dqpsk.theta > PI {
                pi4dqpsk.theta -= 2.0 * PI;
            } else if pi4dqpsk.theta < -PI {
                pi4dqpsk.theta += 2.0 * PI;
            }

            let y = Complex32::from_polar(1.0, pi4dqpsk.theta);
            Ok(y)
        } else {
            Err(Error::Internal("invalid modem data".to_string()))
        }
    }

    pub(super) fn demodulate_pi4dqpsk(&mut self, x: Complex32) -> Result<u32> {
        if let ModemData::Pi4Dqpsk(pi4dqpsk) = self.data.as_mut().unwrap() {
            let theta = x.arg();
            let mut d_theta = theta - pi4dqpsk.theta;
            while d_theta > PI {
                d_theta -= 2.0 * PI;
            }
            while d_theta < -PI {
                d_theta += 2.0 * PI;
            }

            let sym_out = match d_theta {
                d_theta if d_theta > 0.5 * PI => 1,
                d_theta if d_theta > 0.0 => 0,
                d_theta if d_theta < -0.5 * PI => 3,
                _ => 2,
            };

            let d_theta_ideal = match sym_out {
                0 => 0.25 * PI,
                1 => 0.75 * PI,
                2 => -0.25 * PI,
                3 => -0.75 * PI,
                _ => return Err(Error::Internal("invalid output symbol".to_string())),
            };

            self.x_hat = Complex32::from_polar(1.0, pi4dqpsk.theta + d_theta_ideal);
            self.r = x;
            pi4dqpsk.theta = theta;
            Ok(sym_out)
        } else {
            Err(Error::Internal("invalid modem data".to_string()))
        }
    }

    fn demodulate_soft_pi4dqpsk(&mut self, x: Complex32, soft_bits: &mut [u8]) -> Result<u32> {
        let s = self.demodulate_pi4dqpsk(x)?;
        soft_bits[0] = if s & 2 == 0 { 0 } else { 255 };
        soft_bits[1] = if s & 1 == 0 { 0 } else { 255 };
        Ok(s)
    }
}

