import click
from PIL import Image

@click.command()
@click.argument('scale', type=click.INT)
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('output_file', type=click.Path())
def convert_geotiff_to_png(scale, input_file, output_file):
    """Convert a GeoTIFF file to a PNG file."""
    with Image.open(input_file) as img:
        width, height = img.size
        
        rw_width = width * scale
        rw_height = height * scale
        
        _max = max(rw_width, rw_height)
                
        out = Image.new('RGBA', img.size)
        for x in range(width):
            for y in range(height):
                tif_pixel = img.getpixel((x, y))
                out_pixel = int(tif_pixel / _max * 255)
                
                out.putpixel((x,y), (out_pixel, out_pixel, out_pixel, 255))
        
        out.save(output_file, 'PNG')
    click.echo(f"Converted {input_file} to {output_file}")

if __name__ == '__main__':
    convert_geotiff_to_png()