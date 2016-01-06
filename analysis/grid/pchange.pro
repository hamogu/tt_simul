pro pchange,	p_star,z_star, p_ram,shock_extend,$   ;input
  	z_ram,z_end,p_end   ;output

index=sort(p_ram)
z_ram=spline(p_star,z_star,p_ram[index])
z_ram=z_ram[sort(index)]

z_end=z_ram+shock_extend

index=sort(z_end)
p_end=spline(z_star,p_star,z_end[index])
p_end=p_end[sort(index)]

p_ratio=p_end/p_ram

end


