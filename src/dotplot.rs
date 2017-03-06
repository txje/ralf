extern crate image;

use std::fs::File;
use std::path::Path;

use self::image::GenericImage;

pub struct DotPlot {
  pub matrix: Vec<Vec<Vec<Vec<u32>>>>
}

pub fn draw_dp(dp:DotPlot) {
  let gap:usize = 10;

  let mut target_pts:usize = 0;
  let target_blocks:usize = dp.matrix.len();
  for i in 0..target_blocks {
    target_pts += dp.matrix[i][0].len();
  }

  let mut query_pts:usize = 0;
  let query_blocks:usize = dp.matrix[0].len();
  for i in 0..query_blocks {
    query_pts += dp.matrix[0][i][0].len();
  }

  let imgx:u32 = ((target_blocks-1) * gap + target_pts) as u32;
  let imgy:u32 = ((query_blocks-1) * gap + query_pts) as u32;

  debug!("Drawing dotplot -- image size: {}w x {}h", imgx, imgy);

  // Create a new ImgBuf with width: imgx and height: imgy
  let mut im:image::DynamicImage = image::DynamicImage::new_rgb8(imgx, imgy);

  let mut x:u32 = 0;
  let mut y:u32 = 0;
  for t in 0..target_blocks {
    for tpos in 0..dp.matrix[t][0].len() {
      y = 0;
      for q in 0..query_blocks {
        for qpos in 0..dp.matrix[t][q][tpos].len() {
          im.put_pixel(x, y, if dp.matrix[t][q][tpos][qpos] > 0 {image::Rgba{data: [255, 255, 255, 1]}} else {image::Rgba{data: [50, 50, 50, 1]}});
          y += 1;
        }
        y += gap as u32;
      }
      x += 1;
    }
    x += gap as u32;
  }

  let ref mut fout = File::create(&Path::new("dp.png")).unwrap();
  let _ = im.save(fout, image::PNG);
}
